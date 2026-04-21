import polars as pl

"""
Convert BOLD or UNITE search results to format expected by LCA Galaxy tool.
This script reads a BOLD result file in TSV format and converts taxonomy strings
from k__kingdom;p__phylum;c__class format to kingdom / phylum / class .... format.

Usage:
    python bold_to_lca_format.py --blast_result_file <input.tsv> --output_file <output.tsv> --source_db [UNITE|BOLD]
"""
import argparse
from pathlib import Path


def prep_lca(blast_result_file, source):
    header = "#Query ID\t#Subject\t#Subject accession\t#Subject Taxonomy ID\t#Identity percentage\t#Coverage\t#evalue\t#bitscore"
    result = pl.read_csv(
        blast_result_file, separator="\t", new_columns=header.split("\t")
    )

    # add source
    result = result.with_columns(pl.lit(source).alias("#Source"))

    # add Taxonomy
    result = result.with_columns(
        pl.col("#Subject")
        .map_elements(
            lambda s: " / ".join(
                [
                    part.split("__", 1)[1]
                    for part in s.split(";")
                    if part.split("__")[0] in ["k", "p", "c", "o", "f", "g", "s"]
                ]
            ),
            return_dtype=pl.Utf8,
        )
        .alias("#Taxonomy")
    )

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert BOLD or UNITE result to to format expected by LCA tool"
    )
    parser.add_argument("--blast_result_file", help="Path to blast result file")
    parser.add_argument("--output_file", help="Path to formatted output")
    parser.add_argument(
        "--source_db",
        choices=["UNITE", "BOLD"],
        help="Which database was blasted [UNITE, or BOLD]",
    )
    args = parser.parse_args()
    result = prep_lca(args.blast_result_file, source=args.source_db)
    result.write_csv(args.output_file, separator="\t")
