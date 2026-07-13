#!/usr/bin/env python
"""Deterministic runtime and peak-memory benchmarks for OrthoSNAP.

The benchmark deliberately keeps data generation outside measured regions.  It
uses repeated taxon names across otherwise single-copy subgroups, which models
a large multi-copy orthogroup and exercises both candidate discovery and
inparalog checks.  No generated data is retained after the process exits.
"""

import argparse
import gc
import hashlib
import json
import os
import statistics
import subprocess
import sys
import tempfile
import time
import tracemalloc
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from orthosnap.helper import InparalogToKeep
from orthosnap.orthosnap import _extract_subgroups, _write_consensus_outputs, execute


def _combine_balanced(clades):
    """Combine clades into a deterministic, approximately balanced tree."""
    current = list(clades)
    while len(current) > 1:
        next_level = []
        for index in range(0, len(current), 2):
            if index + 1 == len(current):
                next_level.append(current[index])
            else:
                next_level.append(
                    Clade(
                        branch_length=0.1,
                        confidence=100,
                        clades=[current[index], current[index + 1]],
                    )
                )
        current = next_level
    return current[0]


def _combine_unbalanced(clades):
    """Combine clades into a deterministic pectinate tree."""
    iterator = iter(clades)
    current = next(iterator)
    for clade in iterator:
        current = Clade(
            branch_length=0.1,
            confidence=100,
            clades=[current, clade],
        )
    return current


def build_dataset(group_count, taxa_per_group, shape):
    """Return a synthetic multi-copy tree, FASTA mapping, and subgroup tips."""
    combine = _combine_balanced if shape == "balanced" else _combine_unbalanced
    fasta_dict = {}
    subgroup_tips = []
    subgroup_clades = []

    for group_index in range(group_count):
        tips = []
        tip_clades = []
        for taxon_index in range(taxa_per_group):
            name = f"taxon_{taxon_index:04d}|gene_{group_index:05d}"
            tips.append(name)
            tip_clades.append(Clade(name=name, branch_length=0.1))
            fasta_dict[name] = SeqRecord(
                Seq("M" + ("A" * (59 + ((group_index + taxon_index) % 7)))),
                id=name,
                description="",
            )
        subgroup_tips.append(tips)
        subgroup_clades.append(combine(tip_clades))

    return Tree(root=combine(subgroup_clades)), fasta_dict, subgroup_tips


def _measure(function, warmups, runs, measure_memory=True):
    for _ in range(warmups):
        function()

    elapsed = []
    signatures = []
    for _ in range(runs):
        gc.collect()
        started = time.perf_counter()
        signatures.append(function())
        elapsed.append(time.perf_counter() - started)

    peaks = []
    if measure_memory:
        for _ in range(min(runs, 3)):
            gc.collect()
            tracemalloc.start()
            signatures.append(function())
            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            peaks.append(peak / (1024 * 1024))

    if len(set(signatures)) != 1:
        raise RuntimeError(f"benchmark result changed between runs: {signatures}")

    result = {
        "median_seconds": statistics.median(elapsed),
        "min_seconds": min(elapsed),
        "max_seconds": max(elapsed),
        "result_signature": signatures[0],
    }
    if peaks:
        result.update(
            {
                "median_peak_mib": statistics.median(peaks),
                "max_peak_mib": max(peaks),
            }
        )
    return result


def benchmark_core(args):
    tree, fasta_dict, _ = build_dataset(
        args.groups, args.taxa_per_group, args.shape
    )

    def run():
        with open(os.devnull, "w") as devnull:
            with redirect_stdout(devnull), redirect_stderr(devnull):
                result = _extract_subgroups(
                    tree=tree,
                    fasta="synthetic.fasta",
                    fasta_dict=fasta_dict,
                    support=80,
                    occupancy=args.taxa_per_group,
                    snap_trees=False,
                    inparalog_to_keep=InparalogToKeep.longest_seq_len,
                    output_path="",
                    report_inparalog_handling=False,
                    delimiter="|",
                    write_outputs=False,
                )
        groups = ";".join(
            ",".join(record["tips"]) for record in result["subgroup_records"]
        )
        digest = hashlib.sha256(groups.encode("utf-8")).hexdigest()[:16]
        return result["subgroup_counter"], digest

    return _measure(run, args.warmups, args.runs, measure_memory=not args.skip_memory)


def _write_reference_inputs(directory, tree, fasta_dict):
    tree_path = directory / "synthetic.treefile"
    fasta_path = directory / "synthetic.fasta"
    Phylo.write(tree, tree_path, "newick")
    with open(fasta_path, "w") as handle:
        for record in fasta_dict.values():
            handle.write(f">{record.id}\n{record.seq}\n")
    return tree_path, fasta_path


def benchmark_consensus(args):
    tree, fasta_dict, subgroup_tips = build_dataset(
        args.groups, args.taxa_per_group, args.shape
    )
    temporary_directory = tempfile.TemporaryDirectory()
    root = Path(temporary_directory.name)
    tree_path, fasta_path = _write_reference_inputs(root, tree, fasta_dict)
    support_counts = Counter({frozenset(tips): 10 for tips in subgroup_tips})
    run_number = 0

    def run():
        nonlocal run_number
        output_path = root / f"run_{run_number}"
        run_number += 1
        output_path.mkdir()
        with open(os.devnull, "w") as devnull:
            with redirect_stdout(devnull), redirect_stderr(devnull):
                _, emitted = _write_consensus_outputs(
                    fasta=str(fasta_path),
                    fasta_dict=fasta_dict,
                    output_path=f"{output_path}/",
                    delimiter="|",
                    support_counts=support_counts,
                    num_trees=10,
                    min_frequency=0.5,
                    consensus_trees=True,
                    reference_tree_path=str(tree_path),
                    rooted=True,
                )
        tree_outputs = len(list(output_path.glob("*.tre")))
        return emitted, tree_outputs

    try:
        return _measure(run, args.warmups, args.runs, measure_memory=not args.skip_memory)
    finally:
        temporary_directory.cleanup()


def benchmark_bootstrap(args):
    tree, fasta_dict, _ = build_dataset(
        args.groups, args.taxa_per_group, args.shape
    )
    temporary_directory = tempfile.TemporaryDirectory()
    root = Path(temporary_directory.name)
    tree_path, fasta_path = _write_reference_inputs(root, tree, fasta_dict)
    bootstrap_path = root / "bootstrap_trees.txt"
    bootstrap_path.write_text(
        "".join(f"{tree_path}\n" for _ in range(args.bootstrap_count))
    )
    run_number = 0

    def run():
        nonlocal run_number
        output_path = root / f"bootstrap_run_{run_number}"
        run_number += 1
        with open(os.devnull, "w") as devnull:
            with redirect_stdout(devnull), redirect_stderr(devnull):
                result = execute(
                    tree=str(tree_path),
                    fasta=str(fasta_path),
                    support=80,
                    occupancy=args.taxa_per_group,
                    rooted=True,
                    snap_trees=False,
                    inparalog_to_keep=InparalogToKeep.longest_seq_len,
                    report_inparalog_handling=False,
                    output_path=str(output_path),
                    delimiter="|",
                    bootstrap_trees=str(bootstrap_path),
                    consensus_trees=False,
                )
        return result["subgroup_counter"]

    try:
        return _measure(run, args.warmups, args.runs, measure_memory=not args.skip_memory)
    finally:
        temporary_directory.cleanup()


def benchmark_startup(args):
    command = [sys.executable, "-m", "orthosnap", "--help"]

    def run():
        completed = subprocess.run(
            command,
            cwd=Path(__file__).resolve().parents[1],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        return completed.returncode

    return _measure(run, args.warmups, args.runs, measure_memory=False)


def create_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("core", "bootstrap", "consensus", "startup", "all"),
        default="all",
    )
    parser.add_argument("--shape", choices=("balanced", "unbalanced"), default="unbalanced")
    parser.add_argument("--groups", type=int, default=80)
    parser.add_argument("--taxa-per-group", type=int, default=8)
    parser.add_argument("--bootstrap-count", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--skip-memory", action="store_true")
    return parser


def main(argv=None):
    args = create_parser().parse_args(argv)
    benchmarks = {
        "core": benchmark_core,
        "bootstrap": benchmark_bootstrap,
        "consensus": benchmark_consensus,
        "startup": benchmark_startup,
    }
    selected = benchmarks if args.mode == "all" else {args.mode: benchmarks[args.mode]}
    result = {
        "parameters": {
            "groups": args.groups,
            "taxa_per_group": args.taxa_per_group,
            "tips": args.groups * args.taxa_per_group,
            "bootstrap_count": args.bootstrap_count,
            "shape": args.shape,
            "warmups": args.warmups,
            "runs": args.runs,
            "python": sys.version.split()[0],
        },
        "benchmarks": {name: function(args) for name, function in selected.items()},
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
