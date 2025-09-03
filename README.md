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

## GUI (PyQt6)

Install GUI extras and launch the app:

```
pip install -e .[gui]
condrec-gui
# or: python -m condrec.gui
```

UI layout: top pane shows responses; bottom row has an input box for reaction SMILES, a JSON toggle, and buttons for Open, Analyze, Clear.

## Examples

Run built-in examples (after install):

```
python -m condrec analyze --file examples/sn2.txt
python -m condrec analyze --file examples/esterification.txt
python -m condrec analyze --file examples/with_agents.txt
python -m condrec analyze --file examples/invalid.txt
```

You can also pipe input or request JSON:

```
type examples\sn2.txt | python -m condrec analyze --json
# or on Unix:
# cat examples/sn2.txt | python -m condrec analyze --json
```

## Tests

Install dev deps and run pytest:

```
pip install -r requirements-dev.txt
pytest -q
```
