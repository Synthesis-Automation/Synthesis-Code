import os
import sys

import pytest


rdkit = pytest.importorskip("rdkit", reason="RDKit not installed; install via conda-forge")


from condrec.analysis import analyze_reaction_smiles


def test_three_part_reaction_parsing_and_counts():
    rxn = "CCO.CCO>C(N)C=O>CC=O.CCO"
    result = analyze_reaction_smiles(rxn)

    assert isinstance(result, dict)
    assert len(result["reactants"]) == 2
    assert len(result["agents"]) == 1
    assert len(result["products"]) == 2

    # Check all molecules parsed
    assert all(mol.get("ok") for mol in result["reactants"]) is True
    assert all(mol.get("ok") for mol in result["agents"]) is True
    assert all(mol.get("ok") for mol in result["products"]) is True

    # Element count dicts present
    assert isinstance(result["element_counts"]["reactants"], dict)
    assert isinstance(result["element_counts"]["products"], dict)


def test_two_part_normalizes_agents():
    rxn = "CCO>CC=O"  # agents omitted
    result = analyze_reaction_smiles(rxn)
    assert len(result["agents"]) == 0
    assert len(result["reactants"]) == 1
    assert len(result["products"]) == 1


def test_invalid_reaction_raises():
    with pytest.raises(ValueError):
        analyze_reaction_smiles("CCO")  # missing '>' separators


def test_invalid_smiles_reports_ok_false():
    rxn = "C1CCCC1C(>agents)>O"  # malformed reactant; also malformed split, fix to 3 parts
    # Use an intentionally bad SMILES in one position, keep split correct
    rxn = "notasmiles>>O"  # two-part input with invalid reactant
    result = analyze_reaction_smiles(rxn)
    assert len(result["reactants"]) == 1
    assert result["reactants"][0]["ok"] is False

