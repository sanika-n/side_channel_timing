#!/usr/bin/env python3
"""Collect repeated Activity 4 samples and analyze them with the RTLF bootstrap rule."""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import subprocess
import tempfile
import time
import zipfile
import zlib
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
ARTIFACT_ZIP = ROOT / "data" / "f_source_files_artificial_side_channels.zip"
DEFAULT_OUTPUT = ROOT / "activity4_results"
ZENODO_URL = (
    "https://zenodo.org/records/10817685/files/"
    "f_source_files_artificial_side_channels.zip?download=1"
)
CHANNELS = ("shift-1%", "same-xy")
PROBABILITIES = np.arange(1, 10, dtype=float) / 10.0


def get_sources(channel: str) -> tuple[str, str]:
    """Read one server/client pair from the small author-provided source archive."""
    if not ARTIFACT_ZIP.exists():
        raise FileNotFoundError(
            f"Missing {ARTIFACT_ZIP}. Download the 15.7 KB source archive from {ZENODO_URL}."
        )
    with zipfile.ZipFile(ARTIFACT_ZIP) as archive:
        names = archive.namelist()
        folder = f"source-files-{channel}/"
        server_name = next(
            name for name in names if name.endswith(folder + f"server_{channel}.c")
        )
        client_name = next(
            name for name in names if name.endswith(folder + f"client_{channel}.c")
        )
        return (
            archive.read(server_name).decode("utf-8"),
            archive.read(client_name).decode("utf-8"),
        )


def make_bounded_client(source: str) -> str:
    """Add a total-observation argument to the authors' fixed-size client."""
    replacements = {
        "int main(int argc, char const *argv[]) {": (
            "int main(int argc, char const *argv[]) {\n"
            "    int total_samples = (argc > 1) ? atoi(argv[1]) : 20000;\n"
            "    if (total_samples < 2 || total_samples > 1000000 || total_samples % 2 != 0) {\n"
            '        fprintf(stderr, "sample count must be even and between 2 and 1000000\\n");\n'
            "        return 2;\n"
            "    }"
        ),
        "char array[1000000];": "char *array = malloc((size_t)total_samples);\n    if (!array) return 2;",
        "fill_array(array, 1000000);": "fill_array(array, (size_t)total_samples);",
        "for(int i= 0; i<1000000;i++)": "for(int i = 0; i < total_samples; i++)",
    }
    for old, new in replacements.items():
        if old not in source:
            raise ValueError(f"Expected artifact client code was not found: {old!r}")
        source = source.replace(old, new, 1)
    return source


def collect(
    server_bin: Path,
    client_bin: Path,
    count_per_class: int,
    out_csv: Path,
) -> None:
    """Collect balanced BASELINE/MODIFIED round-trip cycle measurements."""
    total = count_per_class * 2
    server_log = out_csv.with_suffix(".server.log")
    with server_log.open("wb") as log:
        server = subprocess.Popen([str(server_bin)], stdout=log, stderr=log)
        try:
            # The source client waits five seconds before connecting, allowing bind().
            time.sleep(0.25)
            if server.poll() is not None:
                raise RuntimeError(
                    f"Artifact server exited early (status {server.returncode}); "
                    f"see {server_log}. Port 8080 may already be in use."
                )
            result = subprocess.run(
                [str(client_bin), str(total)],
                check=True,
                capture_output=True,
                text=True,
                timeout=600,
            )
        finally:
            server.terminate()
            try:
                server.wait(timeout=2)
            except subprocess.TimeoutExpired:
                server.kill()
                server.wait()

    rows = list(csv.DictReader(result.stdout.splitlines()))
    if len(rows) != total:
        raise RuntimeError(f"Expected {total} timing rows, received {len(rows)}.")
    counts = {"BASELINE": 0, "MODIFIED": 0}
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["class", "cycles"])
        writer.writeheader()
        for row in rows:
            label = row["V1"]
            cycles = int(row["V2"])
            if label not in counts or cycles <= 0:
                continue
            writer.writerow({"class": label, "cycles": cycles})
            counts[label] += 1
    if counts != {"BASELINE": count_per_class, "MODIFIED": count_per_class}:
        raise RuntimeError(f"Expected balanced classes; collected {counts}.")


def _type2_quantile(sorted_values: np.ndarray, probability: float) -> float:
    """R quantile(type=2), used by the RTLF implementation."""
    n = len(sorted_values)
    h = n * probability
    nearest = round(h)
    if abs(h - nearest) < 1e-9:
        k = int(nearest)
        if k <= 0:
            return float(sorted_values[0])
        if k >= n:
            return float(sorted_values[-1])
        return float((sorted_values[k - 1] + sorted_values[k]) / 2.0)
    return float(sorted_values[math.ceil(h) - 1])


def _bootstrap_quantiles(
    sorted_data: np.ndarray,
    sample_size: int,
    probability: float,
    repetitions: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Simulate RTLF's bootstrap quantiles via order-statistic beta identities.

    A bootstrap sample of size n from the empirical distribution has its k-th
    order statistic at F_empirical_inverse(U_(k)), where U_(k) is Beta(k,n+1-k).
    For type-2 quantiles at integer n*p, the adjacent order statistics are drawn
    jointly using their conditional spacing. This yields the same bootstrap
    quantile distribution without repeatedly materializing and sorting n values.
    """
    h = sample_size * probability
    nearest = round(h)
    size = len(sorted_data)

    if abs(h - nearest) < 1e-9:
        k = int(nearest)
        u_k = rng.beta(k, sample_size + 1 - k, size=repetitions)
        next_spacing = rng.beta(1, sample_size - k, size=repetitions)
        u_next = u_k + (1.0 - u_k) * next_spacing
        idx_k = np.minimum((u_k * size).astype(np.int64), size - 1)
        idx_next = np.minimum((u_next * size).astype(np.int64), size - 1)
        return (sorted_data[idx_k] + sorted_data[idx_next]) / 2.0

    k = math.ceil(h)
    u_k = rng.beta(k, sample_size + 1 - k, size=repetitions)
    indices = np.minimum((u_k * size).astype(np.int64), size - 1)
    return sorted_data[indices].astype(float)


def rtlf_analysis(
    baseline_values: np.ndarray,
    modified_values: np.ndarray,
    alpha: float,
    bootstrap_repetitions: int,
    seed: int,
) -> dict[str, object]:
    """Apply RTLF's nine-decile, Bonferroni-corrected empirical bootstrap rule."""
    x = np.sort(np.asarray(baseline_values, dtype=np.int64))
    y = np.sort(np.asarray(modified_values, dtype=np.int64))
    n = min(len(x), len(y))
    if n < 100:
        raise ValueError("RTLF requires at least 100 observations per class.")
    x, y = x[:n], y[:n]
    rng = np.random.default_rng(seed)

    observed = []
    thresholds = []
    significant_deciles = []
    alpha_per_decile = alpha / len(PROBABILITIES)
    threshold_probability = 1.0 - alpha_per_decile

    for index, probability in enumerate(PROBABILITIES, start=1):
        observed_difference = abs(
            _type2_quantile(x, float(probability))
            - _type2_quantile(y, float(probability))
        )
        x_boot_1 = _bootstrap_quantiles(x, n, float(probability), bootstrap_repetitions, rng)
        x_boot_2 = _bootstrap_quantiles(x, n, float(probability), bootstrap_repetitions, rng)
        y_boot_1 = _bootstrap_quantiles(y, n, float(probability), bootstrap_repetitions, rng)
        y_boot_2 = _bootstrap_quantiles(y, n, float(probability), bootstrap_repetitions, rng)
        x_noise = np.abs(x_boot_1 - x_boot_2)
        y_noise = np.abs(y_boot_1 - y_boot_2)
        threshold_x = _type2_quantile(np.sort(x_noise), threshold_probability)
        threshold_y = _type2_quantile(np.sort(y_noise), threshold_probability)
        threshold = max(threshold_x, threshold_y)
        is_significant = observed_difference > threshold
        observed.append(
            {
                "decile": index * 10,
                "observed_difference_cycles": observed_difference,
                "threshold_cycles": threshold,
                "threshold_baseline_cycles": threshold_x,
                "threshold_modified_cycles": threshold_y,
                "significant": is_significant,
            }
        )
        thresholds.append(threshold)
        if is_significant:
            significant_deciles.append(index * 10)

    return {
        "n_per_class": n,
        "alpha": alpha,
        "alpha_per_decile": alpha_per_decile,
        "bootstrap_repetitions": bootstrap_repetitions,
        "bootstrap_seed": seed,
        "difference_detected": bool(significant_deciles),
        "significant_deciles": significant_deciles,
        "baseline_median_cycles": _type2_quantile(x, 0.5),
        "modified_median_cycles": _type2_quantile(y, 0.5),
        "decile_results": observed,
    }


def read_groups(path: Path) -> tuple[np.ndarray, np.ndarray]:
    groups: dict[str, list[int]] = {"BASELINE": [], "MODIFIED": []}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            groups[row["class"]].append(int(row["cycles"]))
    return np.asarray(groups["BASELINE"]), np.asarray(groups["MODIFIED"])


def compile_channel(channel: str, work: Path) -> tuple[Path, Path]:
    server_source, client_source = get_sources(channel)
    server_c = work / f"server-{channel.replace('%', 'pct')}.c"
    client_c = work / "client.c"
    server_bin = work / f"server-{channel.replace('%', 'pct')}"
    client_bin = work / "client"
    server_c.write_text(server_source)
    client_c.write_text(make_bounded_client(client_source))
    subprocess.run(["gcc", "-O2", str(server_c), "-o", str(server_bin)], check=True)
    subprocess.run(["gcc", "-O2", str(client_c), "-o", str(client_bin)], check=True)
    return server_bin, client_bin


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sizes",
        type=int,
        nargs="+",
        default=[30000, 200000],
        help="observations per class for each dataset (default: 30000 200000)",
    )
    parser.add_argument("--repeats", type=int, default=3, help="independent collection runs per size")
    parser.add_argument("--bootstrap", type=int, default=10000, help="RTLF bootstrap repetitions")
    parser.add_argument("--alpha", type=float, default=0.09, help="RTLF Type-I error target")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if platform.machine().lower() not in {"x86_64", "amd64"}:
        parser.error("the artifact client uses x86 RDTSC and requires an x86_64 machine")
    if not args.sizes or any(n < 100 for n in args.sizes):
        parser.error("each sample size must be at least 100 observations per class")
    if any(n * 2 > 1_000_000 for n in args.sizes):
        parser.error("the artifact client is bounded to at most 1,000,000 total observations")
    if args.repeats < 2:
        parser.error("use at least two repeats to compare run-to-run variation")
    if args.bootstrap < 1000:
        parser.error("use at least 1,000 bootstrap replicates")
    if not 0 < args.alpha < 1:
        parser.error("alpha must be between 0 and 1")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    run_results: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="activity4-") as tmp:
        work = Path(tmp)
        for channel_index, channel in enumerate(CHANNELS):
            server_bin, client_bin = compile_channel(channel, work)
            channel_slug = channel.replace("%", "-percent")
            for size_index, count in enumerate(args.sizes):
                for repeat in range(1, args.repeats + 1):
                    stem = f"{channel_slug}-n{count}-r{repeat:02d}"
                    csv_path = args.output_dir / f"{stem}.csv"
                    collect(server_bin, client_bin, count, csv_path)
                    baseline, modified = read_groups(csv_path)
                    seed = zlib.crc32(stem.encode("utf-8"))
                    result = rtlf_analysis(
                        baseline,
                        modified,
                        alpha=args.alpha,
                        bootstrap_repetitions=args.bootstrap,
                        seed=seed,
                    )
                    result.update(
                        {
                            "channel": channel,
                            "repeat": repeat,
                            "dataset_file": csv_path.name,
                        }
                    )
                    (args.output_dir / f"{stem}.rtlf.json").write_text(
                        json.dumps(result, indent=2) + "\n"
                    )
                    run_results.append(result)
                    decision = "DIFFERENCE" if result["difference_detected"] else "NO DIFFERENCE"
                    print(
                        f"{channel:8s} n={count:>6,}/class repeat={repeat}: "
                        f"{decision}; significant deciles={result['significant_deciles']}",
                        flush=True,
                    )

    grouped: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    for result in run_results:
        grouped[(str(result["channel"]), int(result["n_per_class"]))].append(result)
    aggregate = []
    for (channel, count), results in grouped.items():
        detections = sum(bool(result["difference_detected"]) for result in results)
        aggregate.append(
            {
                "channel": channel,
                "n_per_class": count,
                "repeats": len(results),
                "difference_detected_runs": detections,
                "no_difference_runs": len(results) - detections,
                "run_results": [
                    {
                        "repeat": result["repeat"],
                        "difference_detected": result["difference_detected"],
                        "significant_deciles": result["significant_deciles"],
                        "baseline_median_cycles": result["baseline_median_cycles"],
                        "modified_median_cycles": result["modified_median_cycles"],
                        "dataset_file": result["dataset_file"],
                    }
                    for result in results
                ],
            }
        )

    report = {
        "experiment": "Activity 4: sample-size effect with repeated artificial-channel measurements",
        "hypothesis": (
            "At the same channel and analysis settings, larger samples increase the "
            "repeatability of detecting the 1% timing shift; same-xy should remain a no-leakage control."
        ),
        "null_hypothesis": "The BASELINE and MODIFIED timing distributions are equal.",
        "alternative_hypothesis": "At least one timing decile differs beyond the bootstrap noise threshold.",
        "dataset_source": "Authors' f_source_files_artificial_side_channels.zip (shift-1% and same-xy source pairs)",
        "analysis": {
            "method": "RTLF nine-decile empirical-bootstrap decision rule, Python implementation",
            "bootstrap_repetitions": args.bootstrap,
            "alpha": args.alpha,
            "quantile_type": 2,
            "threshold": "per-decile 1-alpha/9 quantile of within-class bootstrap differences; max across classes",
            "optimization": "Beta order-statistic identity avoids allocating full bootstrap samples",
        },
        "sample_size_means": "observations per class, with BASELINE and MODIFIED balanced",
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "gcc": subprocess.run(["gcc", "--version"], capture_output=True, text=True, check=True)
            .stdout.splitlines()[0],
        },
        "aggregate_results": aggregate,
        "individual_results": run_results,
    }
    report_path = args.output_dir / "summary_shift1_samexy_rtlf.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nSaved repeated timing datasets and RTLF results under {args.output_dir}")
    print(f"Summary: {report_path}")


if __name__ == "__main__":
    main()
