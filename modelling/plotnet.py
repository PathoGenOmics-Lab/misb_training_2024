from typing import Callable, Tuple, Mapping
from cobra import Model
import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms import bipartite


PLOT_OPTIONS = {
    "edgecolors": "tab:gray",
    "node_size": 300,
    "alpha": 0.9
}
ZERO = 1e-9
EDGE_WIDTH_MAX = 15
EDGE_WIDTH_MIN = 5
EDGE_STYLE = "arc3,rad=0.2"


class NetPlotter:
    
    REACTION = 0
    METABOLITE = 1
    
    def __init__(self, model: Model):
        self.model = model
        self.g = nx.DiGraph()
        for reaction in self.model.reactions:
            self.g.add_node(reaction.id, bipartite=self.REACTION)
            for reactant in reaction.reactants:
                self.g.add_node(reactant.id, bipartite=self.METABOLITE)
                self.g.add_edge(reactant.id, reaction.id)
            for product in reaction.products:
                self.g.add_node(product.id, bipartite=self.METABOLITE)
                self.g.add_edge(reaction.id, product.id)
    
    def _select_all_nodes(self, omit_metabolites: set) -> Tuple[set, set]:
        reaction_nodes = {n for n, d in self.g.nodes(data=True) if d["bipartite"] == self.REACTION}
        metabolite_nodes = set(self.g) - reaction_nodes - omit_metabolites
        return metabolite_nodes, reaction_nodes
    
    def _select_nodes_with_reactions(self, reaction_ids: set, omit_metabolites: set) -> Tuple[set, set]:
        reaction_nodes = {n for n, d in self.g.nodes(data=True) if (d["bipartite"] == self.REACTION) and (n in reaction_ids)}
        metabolite_nodes = set()
        for reaction in (reaction for reaction in self.model.reactions if reaction.id in reaction_nodes):
            for metabolite in reaction.metabolites:
                metabolite_nodes.add(metabolite.id)
        return metabolite_nodes - omit_metabolites, reaction_nodes

    def plot(self, reaction_ids: set = set(), omit_metabolites: set = set(), layout_func: Callable[[nx.Graph], Mapping] = None) -> None:
        # Get types of nodes
        if len(reaction_ids) == 0:
            metabolite_nodes, reaction_nodes = self._select_all_nodes(omit_metabolites)
            g = self.g
        else:
            metabolite_nodes, reaction_nodes = self._select_nodes_with_reactions(reaction_ids, omit_metabolites)
            g = self.g.subgraph(list((reaction_nodes | metabolite_nodes) - omit_metabolites))
        # Set layout
        if layout_func is not None:
            position = layout_func(g)
        else:
            position = nx.kamada_kawai_layout(g)
        # Plot nodes
        nx.draw_networkx_nodes(g, position, nodelist=list(reaction_nodes), node_color="tab:blue", **PLOT_OPTIONS)
        nx.draw_networkx_nodes(g, position, nodelist=list(metabolite_nodes), node_color="tab:orange", **PLOT_OPTIONS)
        nx.draw_networkx_labels(g, position, font_size=10, font_color="black")
        # Plot solid edges
        nx.draw_networkx_edges(g, position, width=1.0, alpha=1.0, connectionstyle=EDGE_STYLE)
        # Plot flux over edges
        solution = self.model.optimize()
        solution_filtered = solution[list(reaction_nodes)]
        nonzero_fluxes = solution_filtered[solution_filtered.abs() > ZERO]
        if nonzero_fluxes.max() != nonzero_fluxes.min():
            scaled_fluxes = EDGE_WIDTH_MIN + (EDGE_WIDTH_MAX - EDGE_WIDTH_MIN) * (nonzero_fluxes - nonzero_fluxes.min()) / (nonzero_fluxes.max() - nonzero_fluxes.min())
        else:
            scaled_fluxes = (nonzero_fluxes - nonzero_fluxes) + EDGE_WIDTH_MAX
        for reaction_id, width in scaled_fluxes.items():
            if reaction_id in reaction_nodes:
                nx.draw_networkx_edges(
                    g,
                    position,
                    edgelist=set(g.out_edges(reaction_id)) | set(g.in_edges(reaction_id)),
                    width=width,
                    alpha=0.5,
                    edge_color="tab:blue",
                    connectionstyle=EDGE_STYLE,
                    arrowstyle="-"
                )
        # Output plot
        plt.tight_layout()
        plt.axis("off")
        plt.show()
