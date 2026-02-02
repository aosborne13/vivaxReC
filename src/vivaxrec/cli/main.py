import argparse

def main():
    parser = argparse.ArgumentParser(
        prog="vivaxrec",
        description="vivaxReC: Vivax Recurrence Consensus"
    )

    subparsers = parser.add_subparsers(dest="command")

    # Placeholder commands
    subparsers.add_parser("build-pairs", help="Build recurrence pair tables")
    subparsers.add_parser("run-aster", help="Run Aster recurrence inference")
    subparsers.add_parser("run-dcifer", help="Run Dcifer relatedness estimation")
    subparsers.add_parser("run-pv3rs", help="Run Pv3Rs classification")
    subparsers.add_parser("consensus", help="Apply vivaxReC consensus logic")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
    else:
        print(f"[vivaxReC] Command '{args.command}' is not implemented yet.")

