# Conditions Recommender (CLI skeleton)

Minimal CLI to accept a reaction SMILES and run a basic RDKit analysis locally.

## Quickstart

1) Create and activate a Python environment (conda recommended for RDKit):

```
conda create -n condrec python=3.10 -y
conda activate condrec
conda install -c conda-forge rdkit -y
```

Alternatively (pip, may be platform-dependent):

```
pip install rdkit-pypi
```

2) Install local package (editable):

```
pip install -e .
```

3) Run a quick analysis:

```
python -m condrec analyze --rxn "CCO.CCO>C(N)C=O>CC=O.CCO"
```

Or pass a file:

```
python -m condrec analyze --file path/to/reaction.txt
```

The CLI parses reaction SMILES as `reactants>agents>products`, analyzes each molecule (formula, MW, rings, aromatic rings, heavy atoms) and prints a summary. If RDKit is missing, the tool will print installation hints.

