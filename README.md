# MISB training 2024

[![PGO badge](https://img.shields.io/badge/PathoGenOmics-Lab-yellow.svg)](https://pathogenomics.github.io/)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/PathoGenOmics-Lab/misb_training_2024/HEAD)

## Running the notebooks locally

### Installing dependencies

An easy way is to create a `conda` environment.
Install it using [`environment.yaml`](/environment.yml),
and then activate it:

```shell
conda env create -f environment.yml
conda activate misb
```

Alternatively, requirements can be installed with `pip`, using
[`requirements.txt`](/requirements.txt):

```shell
pip install -r requirements.txt
```

### Launching Jupyter Lab

Once the dependencies are installed, open the notebooks (`.ipynb` files)
with Jupyter Lab after running the following on the command line:

```shell
jupyter lab
```

## Running the notebooks on a web browser

It's also possible to use [Binder](https://mybinder.org/) to run the
notebooks on your web browser. Click [here](https://mybinder.org/v2/gh/PathoGenOmics-Lab/misb_training_2024/HEAD) or the "launch binder" button at the top!
