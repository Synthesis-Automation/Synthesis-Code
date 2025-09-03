import argparse
import sys
from pathlib import Path

from . import __version__
from .analysis import analyze_reaction_smiles


def _read_reaction_from_file(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8").strip()
    except Exception as exc:
        print(f"Error reading file '{path}': {exc}", file=sys.stderr)
        sys.exit(2)


def _print_analysis(result: dict) -> None:
    # Compact, human-readable summary
    print("=== Reaction Analysis ===")
    print(f"Parsed roles: reactants={len(result['reactants'])}, agents={len(result['agents'])}, products={len(result['products'])}")
    if result.get("atom_map_present"):
        print("Atom maps detected in input.")

    for role in ("reactants", "agents", "products"):
        items = result[role]
        print(f"\n[{role}] count={len(items)}")
        for i, mol in enumerate(items, 1):
            status = "ok" if mol["ok"] else "invalid"
            print(
                f"  {i}. {mol['input_smiles']} -> {status}"
            )
            if mol["ok"]:
                print(
                    f"     formula={mol['formula']} MW={mol['mw']:.2f} heavy={mol['heavy_atoms']} rings={mol['rings']} arom_rings={mol['aromatic_rings']}"
                )
            else:
                print(f"     error={mol['error']}")

    # Simple mass-balance style element accounting (reactants vs products; agents excluded)
    r_counts = result["element_counts"]["reactants"]
    p_counts = result["element_counts"]["products"]
    keys = sorted(set(r_counts) | set(p_counts))
    deltas = {k: p_counts.get(k, 0) - r_counts.get(k, 0) for k in keys}
    print("\n[Element count delta] products - reactants (agents excluded):")
    if not keys:
        print("  (no elements counted)")
    else:
        print(
            "  "
            + ", ".join(
                f"{k}:{deltas[k]:+d}" for k in keys if deltas[k] != 0
            )
            or "  balanced within counted elements"
        )


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="condrec",
        description="Minimal CLI for reaction SMILES analysis (RDKit)",
    )
    parser.add_argument("command", nargs="?", default="analyze", help="Command to run (analyze)")
    parser.add_argument("--version", action="store_true", help="Print version and exit")

    parser.add_argument("--rxn", type=str, help="Reaction SMILES (reactants>agents>products)")
    parser.add_argument("--file", type=str, help="Path to file containing a reaction SMILES")

    args = parser.parse_args(argv)

    if args.version:
        print(__version__)
        return 0

    if args.command not in {"analyze"}:
        print(f"Unknown command: {args.command}", file=sys.stderr)
        return 2

    rxn_smiles = None
    if args.rxn:
        rxn_smiles = args.rxn.strip()
    elif args.file:
        rxn_smiles = _read_reaction_from_file(Path(args.file))

    if not rxn_smiles:
        print("Provide --rxn 'reactants>agents>products' or --file path", file=sys.stderr)
        return 2

    try:
        result = analyze_reaction_smiles(rxn_smiles)
    except ImportError as e:
        print(str(e), file=sys.stderr)
        return 3
    except Exception as e:
        print(f"Failed to analyze reaction: {e}", file=sys.stderr)
        return 1

    _print_analysis(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

