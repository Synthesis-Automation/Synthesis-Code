import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .analysis import analyze_reaction_smiles, format_analysis_summary


def _read_reaction_from_file(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8").strip()
    except Exception as exc:
        print(f"Error reading file '{path}': {exc}", file=sys.stderr)
        sys.exit(2)


def _print_analysis(result: dict) -> None:
    print(format_analysis_summary(result))


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="condrec",
        description="Minimal CLI for reaction SMILES analysis (RDKit)",
    )
    parser.add_argument("command", nargs="?", default="analyze", help="Command to run (analyze)")
    parser.add_argument("--version", action="store_true", help="Print version and exit")

    parser.add_argument("--rxn", type=str, help="Reaction SMILES (reactants>agents>products) or '-' for stdin")
    parser.add_argument("--file", type=str, help="Path to file containing a reaction SMILES or '-' for stdin")
    parser.add_argument("--json", action="store_true", help="Print analysis as JSON")

    args = parser.parse_args(argv)

    if args.version:
        print(__version__)
        return 0

    if args.command not in {"analyze"}:
        print(f"Unknown command: {args.command}", file=sys.stderr)
        return 2

    rxn_smiles = None
    if args.file:
        if args.file == "-":
            rxn_smiles = sys.stdin.read().strip()
        else:
            rxn_smiles = _read_reaction_from_file(Path(args.file))
    elif args.rxn:
        if args.rxn == "-":
            rxn_smiles = sys.stdin.read().strip()
        else:
            rxn_smiles = args.rxn.strip()
    elif not sys.stdin.isatty():
        # Allow piping input without flags
        rxn_smiles = sys.stdin.read().strip()

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

    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        _print_analysis(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
