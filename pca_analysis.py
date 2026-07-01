#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def find_header_line(lines: list[str]) -> int:
    for index, line in enumerate(lines):
        fields = line.split()
        if len(fields) >= 3 and fields[0] == "Name" and fields[1] == "Description":
            return index

    raise ValueError(
        "Could not find a header line starting with 'Name Description'."
    )


def load_count_file(path: str) -> pd.DataFrame:
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]

    header_index = find_header_line(lines)
    data_lines = lines[header_index:]

    if len(data_lines) < 2:
        raise ValueError("File does not contain any count rows after the header.")

    header = data_lines[0].split()
    if len(header) < 3:
        raise ValueError(
            "Header line must contain at least Name, Description, and one sample."
        )

    sample_names = header[2:]
    feature_names = []
    records = []

    for line in data_lines[1:]:
        parts = line.split()
        if len(parts) < 2 + len(sample_names):
            continue

        # The first feature identifier column duplicates the second one.
        feature_name = parts[1]
        counts = [float(value) for value in parts[2 : 2 + len(sample_names)]]

        feature_names.append(feature_name)
        records.append(counts)

    return pd.DataFrame(records, index=feature_names, columns=sample_names)


def run_pca(counts: pd.DataFrame, output_prefix: str = "pca") -> None:
    # samples should be rows for sklearn PCA
    sample_matrix = counts.T

    # optional: remove features with zero variance
    sample_matrix = sample_matrix.loc[:, sample_matrix.var(axis=0) > 0]

    scaler = StandardScaler()
    scaled = scaler.fit_transform(sample_matrix)

    pca = PCA(n_components=2)
    coords = pca.fit_transform(scaled)

    pca_df = pd.DataFrame(
        coords,
        index=sample_matrix.index,
        columns=["PC1", "PC2"],
    )
    pca_df.index.name = "sample"

    explained = pca.explained_variance_ratio_ * 100.0

    plt.figure(figsize=(8, 6))
    plt.scatter(pca_df["PC1"], pca_df["PC2"])

    for sample_name, row in pca_df.iterrows():
        plt.text(row["PC1"], row["PC2"], str(sample_name), fontsize=9)

    plt.xlabel(f"PC1 ({explained[0]:.2f}% variance)")
    plt.ylabel(f"PC2 ({explained[1]:.2f}% variance)")
    plt.title("PCA of count matrix")
    plt.tight_layout()
    plt.savefig(f"{output_prefix}.png", dpi=200)

    pca_df.reset_index().to_csv(f"{output_prefix}_coordinates.tsv", sep="\t", index=False)
    print(f"Saved plot to {output_prefix}.png")
    print(f"Saved PCA coordinates to {output_prefix}_coordinates.tsv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run PCA on a count matrix.")
    parser.add_argument("input_file", help="Path to the count matrix file.")
    parser.add_argument(
        "--output-prefix",
        default="pca",
        help="Prefix for the PCA plot and coordinates output files.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_file = Path(args.input_file).expanduser().resolve()
    counts = load_count_file(str(input_file))
    run_pca(counts, output_prefix=args.output_prefix)


if __name__ == "__main__":
    main()
