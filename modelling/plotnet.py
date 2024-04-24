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
    
    def plot(self, layout_func=None):
        if layout_func is not None:
            position = layout_func(self.g)
        else:
            position = nx.kamada_kawai_layout(self.g)
        # Get types of nodes
        metabolite_nodes = {n for n, d in self.g.nodes(data=True) if d["bipartite"] == self.METABOLITE}
        reaction_nodes = set(self.g) - metabolite_nodes
        # Plot nodes
        nx.draw_networkx_nodes(self.g, position, nodelist=list(reaction_nodes), node_color="tab:blue", **PLOT_OPTIONS)
        nx.draw_networkx_nodes(self.g, position, nodelist=list(metabolite_nodes), node_color="tab:orange", **PLOT_OPTIONS)
        nx.draw_networkx_labels(self.g, position, font_size=10, font_color="white")
        # Plot solid edges
        nx.draw_networkx_edges(self.g, position, width=1.0, alpha=1., connectionstyle=EDGE_STYLE)
        # Plot flux over edges
        solution = self.model.optimize()
        nonzero_fluxes = solution.fluxes[solution.fluxes.abs() > ZERO]
        if nonzero_fluxes.max() != nonzero_fluxes.min():
            scaled_fluxes = EDGE_WIDTH_MIN + (EDGE_WIDTH_MAX - EDGE_WIDTH_MIN) * (nonzero_fluxes - nonzero_fluxes.min()) / (nonzero_fluxes.max() - nonzero_fluxes.min())
        else:
            scaled_fluxes = (nonzero_fluxes - nonzero_fluxes) + nonzero_fluxes.max()
        for reaction_id, width in scaled_fluxes.items():
            # print(f"Adding flux {flux:g} to reaction {reaction_id}")
            nx.draw_networkx_edges(
                self.g,
                position,
                edgelist=set(self.g.out_edges(reaction_id)) | set(self.g.in_edges(reaction_id)),
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
