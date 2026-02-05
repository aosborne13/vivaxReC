import argparse

from vivaxrec.commands.prepare_input_data import run_prepare_input_data_single_dataset


def main():
    parser = argparse.ArgumentParser(
        prog="vivaxrec",
        description="vivaxReC: Vivax Recurrence Consensus",
    )
    subparsers = parser.add_subparsers(dest="command")

    # ---- prepare-input-data (single dataset) ----
    p = subparsers.add_parser(
        "prepare-input-data",
        help="Format genotype + metadata into Aster-ready inputs for ONE dataset.",
    )
    p.add_argument("--dataset-name", required=True, help="Dataset name label (e.g. 'improv', 'opra', 'siteA_2024').")
    p.add_argument("--genotypes", required=True, help="Genotype TSV (columns: sample_id, locus, allele).")
    p.add_argument("--metadata", required=True, help="Metadata TSV (must include sample_id, day, time_point).")
    p.add_argument("--outdir", required=True, help="Output directory to write formatted inputs.")

    # Optional: tell vivaxReC which metadata column contains dataset/trial (so we can filter)
    p.add_argument(
        "--metadata-dataset-col",
        default=None,
        help="Optional metadata column containing dataset/trial labels (e.g. 'trial' or 'dataset'). "
             "If provided, metadata will be filtered to --dataset-name.",
    )

    # Optional COI source
    p.add_argument(
        "--coi-source",
        default=None,
        help="Optional COI/MOI TSV to merge (if omitted, vivaxReC uses a naive lower bound from genotypes).",
    )
    p.add_argument(
        "--coi-cols",
        nargs="*",
        default=None,
        help="Two column names in the COI file: SAMPLE_COL COI_COL (e.g. 'sample_id post_effective_coi_mean').",
    )

    # Optional parsing controls (keep defaults)
    p.add_argument("--sample-id-col", default="sample_id", help="Genotype+metadata sample ID column name.")
    p.add_argument("--locus-col", default="locus", help="Genotype locus column name.")
    p.add_argument("--allele-col", default="allele", help="Genotype allele column name.")
    p.add_argument("--day-col", default="day", help="Metadata day column name.")
    p.add_argument("--time-point-col", default="time_point", help="Metadata time_point column name.")
    p.add_argument("--patient-id-col", default="patient_id", help="Metadata patient_id column name (if present).")
    p.add_argument("--patient-split-char", default="_", help="If patient_id is missing, split sample_id on this char.")

    # Keep stubs for later
    subparsers.add_parser("calculate-priors", help="(stub) Calculate priors/background parameters")
    subparsers.add_parser("run-aster-dcifer", help="(stub) Run Aster + Dcifer")
    subparsers.add_parser("run-consensus-pipeline", help="(stub) Run full vivaxReC pipeline")
    subparsers.add_parser("plot-outputs", help="(stub) Plot outputs")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    if args.command == "prepare-input-data":
        run_prepare_input_data_single_dataset(
            dataset_name=args.dataset_name,
            genotypes_file=args.genotypes,
            metadata_file=args.metadata,
            outdir=args.outdir,
            coi_source_file=args.coi_source,
            coi_cols=args.coi_cols,
            metadata_dataset_col=args.metadata_dataset_col,
            sample_id_col=args.sample_id_col,
            locus_col=args.locus_col,
            allele_col=args.allele_col,
            day_col=args.day_col,
            time_point_col=args.time_point_col,
            patient_id_col=args.patient_id_col,
            split_char=args.patient_split_char,
        )
        return

    print(f"[vivaxReC] Command '{args.command}' is not implemented yet.")


if __name__ == "__main__":
    main()


