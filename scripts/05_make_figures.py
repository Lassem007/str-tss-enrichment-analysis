#!/usr/bin/env python3
from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RESULTS = Path("results")
FIGURES = Path("figures")
FIGURES.mkdir(exist_ok=True)


def read_metric_table(path: Path) -> dict:
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df.iloc[:, 0], df.iloc[:, 1]))


def fmt_int(x):
    return f"{int(round(float(x))):,}"


def figure_inside_outside():
    path = RESULTS / "str_overlap_summary.tsv"
    if not path.exists():
        print(f"Skipping inside/outside figure: missing {path}")
        return
    m = read_metric_table(path)
    inside = int(float(m["str_loci_inside_pm2kb_tss_windows"]))
    outside = int(float(m["str_loci_outside_pm2kb_tss_windows"]))
    total = inside + outside
    labels = ["Inside ±2 kb\nTSS windows", "Outside ±2 kb\nTSS windows"]
    counts = [inside, outside]
    plt.figure(figsize=(7, 5))
    bars = plt.bar(labels, counts)
    plt.ylabel("Number of STR loci")
    plt.title("STR loci inside versus outside ±2 kb TSS windows")
    for bar, count in zip(bars, counts):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height(), fmt_int(count),
                 ha="center", va="bottom")
    text = f"Inside = {inside/total*100:.2f}%\nOutside = {outside/total*100:.2f}%"
    plt.text(0.5, max(counts)*0.80, text, ha="center", va="center",
             bbox=dict(facecolor="white", edgecolor="black", alpha=0.8))
    plt.tight_layout()
    plt.savefig(FIGURES / "figure_inside_vs_outside_tss_pm2kb.png", dpi=300)
    plt.close()


def figure_distance_bins(bin_size=250, max_distance=2000):
    path = RESULTS / "str_to_nearest_tss.tsv"
    if not path.exists():
        print(f"Skipping distance bins figure: missing {path}")
        return
    distances = []
    with path.open() as f:
        for line in f:
            try:
                d = int(line.rstrip("\n").split("\t")[-1])
            except ValueError:
                continue
            d = abs(d)
            if d <= max_distance:
                distances.append(d)
    bins = list(range(0, max_distance + bin_size, bin_size))
    labels = [f"{bins[i]}-{bins[i+1]} bp" for i in range(len(bins)-1)]
    counts = []
    for i in range(len(bins)-1):
        left, right = bins[i], bins[i+1]
        if i == 0:
            c = sum(left <= d <= right for d in distances)
        else:
            c = sum(left < d <= right for d in distances)
        counts.append(c)
    pd.DataFrame({"distance_bin": labels, "count": counts}).to_csv(
        RESULTS / "distance_bins_250bp.tsv", sep="\t", index=False
    )
    plt.figure(figsize=(8, 5))
    plt.bar(labels, counts)
    plt.xticks(rotation=45, ha="right")
    plt.xlabel("Distance to nearest transcript-level TSS")
    plt.ylabel("Number of STR-associated output lines")
    plt.title("STR distance to nearest TSS in 250 bp bins")
    plt.tight_layout()
    plt.savefig(FIGURES / "figure_distance_250bp_bins.png", dpi=300)
    plt.close()


def figure_coding_vs_noncoding(bin_size=100, limit=2000):
    coding_bed = RESULTS / "str_coding.bed"
    noncoding_bed = RESULTS / "str_noncoding.bed"
    coding_file = RESULTS / "coding_to_nearest_tss_signed.tsv"
    noncoding_file = RESULTS / "noncoding_to_nearest_tss_signed.tsv"
    needed = [coding_bed, noncoding_bed, coding_file, noncoding_file]
    if not all(p.exists() for p in needed):
        print("Skipping coding/non-coding figure: run scripts/04_coding_noncoding_analysis.sh first")
        return

    def count_lines(path):
        with path.open() as f:
            return sum(1 for _ in f)

    def read_signed_distances(path):
        values = []
        with path.open() as f:
            for line in f:
                try:
                    d = int(line.rstrip("\n").split("\t")[-1])
                except ValueError:
                    continue
                if -limit <= d <= limit:
                    values.append(d)
        return np.array(values)

    n_coding = count_lines(coding_bed)
    n_noncoding = count_lines(noncoding_bed)
    coding_dist = read_signed_distances(coding_file)
    noncoding_dist = read_signed_distances(noncoding_file)

    bins = np.arange(-limit, limit + bin_size, bin_size)
    coding_counts, edges = np.histogram(coding_dist, bins=bins)
    noncoding_counts, _ = np.histogram(noncoding_dist, bins=bins)
    coding_norm = coding_counts / n_coding * 10000
    noncoding_norm = noncoding_counts / n_noncoding * 10000
    x = edges[:-1]

    pd.DataFrame({
        "bin_start": x,
        "bin_end": edges[1:],
        "coding_count": coding_counts,
        "noncoding_count": noncoding_counts,
        "coding_normalized_per_10000": coding_norm,
        "noncoding_normalized_per_10000": noncoding_norm,
    }).to_csv(RESULTS / "coding_vs_noncoding_tss_distance_100bp.tsv", sep="\t", index=False)

    plt.figure(figsize=(9, 5))
    plt.step(x, noncoding_norm, where="post", label="Non-coding STRs")
    plt.step(x, coding_norm, where="post", label="Coding STRs")
    plt.axvline(0, linestyle="--")
    plt.xlabel("Distance to nearest transcript-level TSS (bp)")
    plt.ylabel("Normalized counts per 100 bp per 10,000 STR loci")
    plt.title("Coding versus non-coding STRs relative to nearest TSS")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIGURES / "figure_coding_vs_noncoding_tss_2kb_100bp_normalized.png", dpi=300)
    plt.close()


def figure_observed_vs_expected(n=100):
    summary_path = RESULTS / f"enrichment_summary_{n}.tsv"
    if not summary_path.exists():
        candidates = sorted(RESULTS.glob("enrichment_summary_*.tsv"))
        if candidates:
            summary_path = candidates[0]
        else:
            print("Skipping enrichment figure: missing enrichment_summary_*.tsv")
            return
    m = read_metric_table(summary_path)
    observed = float(m["observed_overlap"])
    expected = float(m["expected_mean"])
    fold = float(m["fold_enrichment"])
    pval = float(m["empirical_p_value"])
    z = float(m["z_score"])
    nperm = int(float(m["number_of_permutations"]))

    labels = ["Random mean\n(chromosome-aware)", "Observed\nSTR loci"]
    values = [expected, observed]
    plt.figure(figsize=(7, 5))
    bars = plt.bar(labels, values)
    plt.ylabel("Number of STRs overlapping ±2 kb TSS windows")
    plt.title("Observed vs expected STR overlap near TSS")
    for bar, value in zip(bars, values):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height(), fmt_int(value),
                 ha="center", va="bottom")
    text = f"Permutations = {nperm}\nFold enrichment = {fold:.3f}\nEmpirical p = {pval:.4f}\nZ-score = {z:.2f}"
    plt.text(0.05, max(values)*0.87, text, ha="left", va="top",
             bbox=dict(facecolor="white", edgecolor="black", alpha=0.8))
    plt.tight_layout()
    plt.savefig(FIGURES / "figure_observed_vs_expected_str_overlap_near_tss.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    figure_inside_outside()
    figure_distance_bins()
    figure_coding_vs_noncoding()
    figure_observed_vs_expected()
    print("Figure generation completed.")
