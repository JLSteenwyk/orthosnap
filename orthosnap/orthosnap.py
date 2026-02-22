#!/usr/bin/env python

import csv
import hashlib
import json
import os
import re
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from Bio import Phylo
from Bio import SeqIO
from tqdm import tqdm

from .args_processing import determine_occupancy_threshold, process_args
from .helper import (
    build_subtree_taxa_cache,
    check_if_single_copy,
    get_all_tips_and_taxa_names,
    handle_multi_copy_subtree,
    handle_single_copy_subtree,
    read_input_files,
)
from .helper import InparalogToKeep
from .parser import create_parser
from .plotter import plot_snap_ogs
from .version import __version__
from .writer import write_output_stats, write_user_args


def _sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_bool(value, default=False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _validate_inputs(tree_path: str, fasta_path: str, delimiter: str):
    errors = []

    try:
        tree = Phylo.read(tree_path, "newick")
    except Exception as exc:  # pragma: no cover - defensive
        return False, {"errors": [f"Failed to parse tree: {exc}"]}

    fasta_records = list(SeqIO.parse(fasta_path, "fasta"))
    if not fasta_records:
        errors.append("Input FASTA contains no sequences.")

    fasta_ids = [record.id for record in fasta_records]
    duplicate_ids = [name for name, count in Counter(fasta_ids).items() if count > 1]
    if duplicate_ids:
        errors.append(
            "Duplicate FASTA IDs detected: " + ", ".join(sorted(duplicate_ids)[:10])
        )

    tree_tips = [tip.name for tip in tree.get_terminals() if tip.name is not None]
    missing_tree_names = ["<unnamed tip>"] if len(tree_tips) != len(tree.get_terminals()) else []
    if missing_tree_names:
        errors.append("Tree contains unnamed tips.")

    fasta_missing_delimiter = [seq_id for seq_id in fasta_ids if delimiter not in seq_id]
    tree_missing_delimiter = [tip for tip in tree_tips if delimiter not in tip]

    if fasta_missing_delimiter:
        errors.append(
            f"Delimiter '{delimiter}' missing from {len(fasta_missing_delimiter)} FASTA header(s)."
        )
    if tree_missing_delimiter:
        errors.append(
            f"Delimiter '{delimiter}' missing from {len(tree_missing_delimiter)} tree tip label(s)."
        )

    fasta_set = set(fasta_ids)
    tree_set = set(tree_tips)
    missing_in_tree = sorted(fasta_set - tree_set)
    missing_in_fasta = sorted(tree_set - fasta_set)

    if missing_in_tree:
        errors.append(
            f"{len(missing_in_tree)} FASTA sequence IDs are missing from tree tips."
        )
    if missing_in_fasta:
        errors.append(
            f"{len(missing_in_fasta)} tree tips are missing from FASTA headers."
        )

    taxa = set()
    for seq_id in fasta_ids:
        taxa.add(seq_id.split(delimiter, 1)[0] if delimiter in seq_id else seq_id)

    summary = {
        "tree_tips": len(tree_tips),
        "fasta_sequences": len(fasta_ids),
        "unique_taxa": len(taxa),
        "missing_in_tree": len(missing_in_tree),
        "missing_in_fasta": len(missing_in_fasta),
        "errors": errors,
    }
    return len(errors) == 0, summary


def _write_structured_outputs(
    fasta: str,
    tree: str,
    output_path: str,
    subgroup_records: list,
    start_time: float,
    end_time: float,
    args_snapshot: dict,
    status: str = "completed",
    extra: dict = None,
):
    fasta_path_stripped = re.sub("^.*/", "", fasta)
    prefix = f"{output_path}{fasta_path_stripped}.orthosnap"
    json_path = f"{prefix}.run.json"
    tsv_path = f"{prefix}.subgroups.tsv"

    record_rows = []
    for record in subgroup_records:
        tips = record.get("tips", [])
        taxa = {tip.split(args_snapshot.get("delimiter", "|"), 1)[0] for tip in tips}
        record_rows.append(
            {
                "subgroup_id": record.get("subgroup_id"),
                "tip_count": len(tips),
                "taxa_count": len(taxa),
                "tips": tips,
            }
        )

    with open(tsv_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["subgroup_id", "tip_count", "taxa_count", "tips"])
        for row in record_rows:
            writer.writerow(
                [
                    row["subgroup_id"],
                    row["tip_count"],
                    row["taxa_count"],
                    ";".join(row["tips"]),
                ]
            )

    payload = {
        "status": status,
        "orthosnap_version": __version__,
        "started_at": datetime.fromtimestamp(start_time, tz=timezone.utc).isoformat(),
        "finished_at": datetime.fromtimestamp(end_time, tz=timezone.utc).isoformat(),
        "execution_seconds": round(end_time - start_time, 6),
        "input": {
            "tree": tree,
            "fasta": fasta,
            "tree_sha256": _sha256(tree),
            "fasta_sha256": _sha256(fasta),
        },
        "arguments": args_snapshot,
        "summary": {
            "subgroup_count": len(subgroup_records),
        },
        "outputs": {
            "subgroups_tsv": tsv_path,
        },
        "subgroups": record_rows,
    }
    if extra:
        payload["extra"] = extra

    with open(json_path, "w") as handle:
        json.dump(payload, handle, indent=2)

    return json_path, tsv_path


def _extract_subgroups(
    tree,
    fasta: str,
    fasta_dict: dict,
    support: float,
    occupancy: float,
    snap_trees: bool,
    inparalog_to_keep: InparalogToKeep,
    output_path: str,
    report_inparalog_handling: bool,
    delimiter: str,
    write_outputs: bool,
):
    taxa, all_tips = get_all_tips_and_taxa_names(tree, delimiter)

    if check_if_single_copy(taxa, all_tips):
        return {
            "single_copy": True,
            "subgroup_counter": 0,
            "subgroup_records": [],
        }

    assigned_tips = set()
    subgroup_counter = 0

    inparalog_handling = dict()
    inparalog_handling_summary = dict()
    subgroup_records = []
    subtree_cache = build_subtree_taxa_cache(tree, delimiter)

    for inter in tqdm(tree.get_nonterminals()[1:]):
        (
            terms,
            terms_set,
            counts_of_taxa_from_terms,
            counts,
        ) = subtree_cache[inter]

        if len(counts_of_taxa_from_terms) >= occupancy:
            if set([1]) == set(counts) and assigned_tips.isdisjoint(terms_set):
                (
                    subgroup_counter,
                    assigned_tips,
                    inparalog_handling,
                    inparalog_handling_summary,
                ) = handle_single_copy_subtree(
                    inter,
                    terms,
                    subgroup_counter,
                    fasta,
                    support,
                    fasta_dict,
                    assigned_tips,
                    snap_trees,
                    output_path,
                    inparalog_handling,
                    inparalog_handling_summary,
                    report_inparalog_handling,
                    subgroup_records,
                    write_outputs,
                )
            elif assigned_tips.isdisjoint(terms_set):
                (
                    subgroup_counter,
                    assigned_tips,
                    inparalog_handling,
                    inparalog_handling_summary,
                ) = handle_multi_copy_subtree(
                    inter,
                    terms,
                    subgroup_counter,
                    fasta,
                    support,
                    fasta_dict,
                    assigned_tips,
                    counts_of_taxa_from_terms,
                    snap_trees,
                    inparalog_to_keep,
                    output_path,
                    inparalog_handling,
                    inparalog_handling_summary,
                    report_inparalog_handling,
                    delimiter,
                    subgroup_records,
                    write_outputs,
                )

    return {
        "single_copy": False,
        "subgroup_counter": subgroup_counter,
        "subgroup_records": subgroup_records,
    }


def _load_bootstrap_trees(bootstrap_tree_file: str):
    trees = []
    with open(bootstrap_tree_file, "r") as handle:
        for line in handle:
            value = line.strip()
            if not value or value.startswith("#"):
                continue
            trees.append(value)
    return trees


def _write_consensus_outputs(
    fasta: str,
    fasta_dict: dict,
    output_path: str,
    delimiter: str,
    support_counts: Counter,
    num_trees: int,
    min_frequency: float,
):
    fasta_path_stripped = re.sub("^.*/", "", fasta)
    tsv_path = f"{output_path}{fasta_path_stripped}.orthosnap.consensus.tsv"

    sorted_items = sorted(
        support_counts.items(), key=lambda item: (-item[1], -len(item[0]), tuple(sorted(item[0])))
    )

    emitted = 0
    with open(tsv_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["consensus_id", "count", "frequency", "tip_count", "taxa_count", "tips"])
        for idx, (tips_set, count) in enumerate(sorted_items):
            frequency = count / num_trees
            if frequency < min_frequency:
                continue
            emitted += 1
            tips = sorted(tips_set)
            taxa_count = len({tip.split(delimiter, 1)[0] for tip in tips})
            consensus_id = f"consensus_{idx}"
            writer.writerow([consensus_id, count, round(frequency, 6), len(tips), taxa_count, ";".join(tips)])

            fasta_out = f"{output_path}{fasta_path_stripped}.orthosnap.{consensus_id}.fa"
            with open(fasta_out, "w") as out_handle:
                for tip in tips:
                    if tip in fasta_dict:
                        SeqIO.write(fasta_dict[tip], out_handle, "fasta")

    return tsv_path, emitted


def execute(
    tree: str,
    fasta: str,
    support: float,
    occupancy: float,
    rooted: bool,
    snap_trees: bool,
    inparalog_to_keep: InparalogToKeep,
    report_inparalog_handling: bool,
    output_path: str,
    delimiter: str,
    plot_snap_ogs_output: bool = False,
    plot_format: str = "png",
    occupancy_mode: str = "legacy",
    occupancy_count: int = None,
    occupancy_fraction: float = None,
    validate_only: bool = False,
    resume: bool = False,
    structured_output: bool = False,
    bootstrap_trees: str = None,
    consensus_min_frequency: float = 0.5,
    total_taxa: int = None,
):
    """
    Master execute Function
    -----------------------
    This function executes the main functions and calls other subfunctions
    """
    if not output_path.endswith("/"):
        output_path = output_path + "/"

    os.makedirs(output_path, exist_ok=True)

    valid, validation_summary = _validate_inputs(tree, fasta, delimiter)
    if not valid:
        print("Input validation failed:")
        for error in validation_summary["errors"]:
            print(f"- {error}")
        if structured_output:
            now = time.time()
            _write_structured_outputs(
                fasta=fasta,
                tree=tree,
                output_path=output_path,
                subgroup_records=[],
                start_time=now,
                end_time=now,
                args_snapshot={
                    "support": support,
                    "occupancy": occupancy,
                    "delimiter": delimiter,
                    "validate_only": validate_only,
                },
                status="validation_failed",
                extra={"validation": validation_summary},
            )
        sys.exit(1)

    print("Input validation summary:")
    print(f"- Tree tips: {validation_summary['tree_tips']}")
    print(f"- FASTA sequences: {validation_summary['fasta_sequences']}")
    print(f"- Unique taxa: {validation_summary['unique_taxa']}")

    if validate_only:
        print("Validation-only mode requested; exiting without SNAP-OG extraction.")
        return {
            "status": "validated",
            "subgroup_counter": 0,
            "subgroup_records": [],
        }

    fasta_path_stripped = re.sub("^.*/", "", fasta)
    run_json_path = f"{output_path}{fasta_path_stripped}.orthosnap.run.json"
    run_marker = Path(run_json_path)

    if resume:
        if run_marker.exists():
            try:
                with open(run_marker, "r") as handle:
                    run_payload = json.load(handle)
                if run_payload.get("status") == "completed":
                    print(f"Resume enabled: existing completed run found at {run_json_path}; skipping.")
                    return {
                        "status": "skipped",
                        "subgroup_counter": run_payload.get("summary", {}).get("subgroup_count", 0),
                        "subgroup_records": run_payload.get("subgroups", []),
                    }
            except (json.JSONDecodeError, OSError):
                pass
        else:
            existing = list(Path(output_path).glob(f"{fasta_path_stripped}.orthosnap.*.fa"))
            if existing:
                print("Resume enabled: existing subgroup FASTA outputs detected; skipping.")
                return {
                    "status": "skipped",
                    "subgroup_counter": len(existing),
                    "subgroup_records": [],
                }

    if report_inparalog_handling:
        inparalog_report_output_name = fasta_path_stripped + ".inparalog_report.txt"
        if os.path.isfile(f"{output_path}{inparalog_report_output_name}"):
            os.remove(f"{output_path}{inparalog_report_output_name}")

    write_user_args(
        tree,
        fasta,
        support,
        occupancy,
        rooted,
        snap_trees,
        inparalog_to_keep,
        report_inparalog_handling,
        output_path,
        delimiter,
        plot_snap_ogs_output,
        plot_format,
    )

    start_time = time.time()

    if bootstrap_trees:
        tree_paths = _load_bootstrap_trees(bootstrap_trees)
        if not tree_paths:
            print("No bootstrap tree paths were provided.")
            sys.exit(1)

        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        support_counts = Counter()
        for tree_path in tree_paths:
            tree_obj, _ = read_input_files(tree_path, fasta, rooted)
            extraction = _extract_subgroups(
                tree=tree_obj,
                fasta=fasta,
                fasta_dict=fasta_dict,
                support=support,
                occupancy=occupancy,
                snap_trees=False,
                inparalog_to_keep=inparalog_to_keep,
                output_path=output_path,
                report_inparalog_handling=False,
                delimiter=delimiter,
                write_outputs=False,
            )
            subgroup_sets = {frozenset(record["tips"]) for record in extraction["subgroup_records"]}
            for subgroup in subgroup_sets:
                support_counts[subgroup] += 1

        consensus_tsv, emitted = _write_consensus_outputs(
            fasta=fasta,
            fasta_dict=fasta_dict,
            output_path=output_path,
            delimiter=delimiter,
            support_counts=support_counts,
            num_trees=len(tree_paths),
            min_frequency=consensus_min_frequency,
        )
        end_time = time.time()

        print(f"Consensus groups meeting frequency >= {consensus_min_frequency}: {emitted}")
        print(f"Consensus summary TSV: {consensus_tsv}")

        if structured_output:
            _write_structured_outputs(
                fasta=fasta,
                tree=tree,
                output_path=output_path,
                subgroup_records=[],
                start_time=start_time,
                end_time=end_time,
                args_snapshot={
                    "support": support,
                    "occupancy": occupancy,
                    "occupancy_mode": occupancy_mode,
                    "occupancy_count": occupancy_count,
                    "occupancy_fraction": occupancy_fraction,
                    "rooted": rooted,
                    "delimiter": delimiter,
                    "bootstrap_trees": bootstrap_trees,
                    "consensus_min_frequency": consensus_min_frequency,
                },
                status="completed",
                extra={
                    "consensus_tsv": consensus_tsv,
                    "consensus_groups_emitted": emitted,
                    "bootstrap_tree_count": len(tree_paths),
                },
            )

        return {
            "status": "completed",
            "subgroup_counter": emitted,
            "subgroup_records": [],
        }

    tree_obj, fasta_dict = read_input_files(tree, fasta, rooted)

    extraction = _extract_subgroups(
        tree=tree_obj,
        fasta=fasta,
        fasta_dict=fasta_dict,
        support=support,
        occupancy=occupancy,
        snap_trees=snap_trees,
        inparalog_to_keep=inparalog_to_keep,
        output_path=output_path,
        report_inparalog_handling=report_inparalog_handling,
        delimiter=delimiter,
        write_outputs=True,
    )

    subgroup_counter = extraction["subgroup_counter"]
    subgroup_records = extraction["subgroup_records"]

    plot_file = None
    if plot_snap_ogs_output and subgroup_counter > 0:
        plot_file = plot_snap_ogs(
            tree=tree_obj,
            subgroup_records=subgroup_records,
            fasta=fasta,
            output_path=output_path,
            plot_format=plot_format,
        )

    write_output_stats(
        fasta,
        subgroup_counter,
        start_time,
        snap_trees,
        output_path,
        plot_file,
    )

    end_time = time.time()

    if structured_output:
        _write_structured_outputs(
            fasta=fasta,
            tree=tree,
            output_path=output_path,
            subgroup_records=subgroup_records,
            start_time=start_time,
            end_time=end_time,
            args_snapshot={
                "support": support,
                "occupancy": occupancy,
                "occupancy_mode": occupancy_mode,
                "occupancy_count": occupancy_count,
                "occupancy_fraction": occupancy_fraction,
                "rooted": rooted,
                "snap_trees": snap_trees,
                "inparalog_to_keep": inparalog_to_keep.value,
                "report_inparalog_handling": report_inparalog_handling,
                "output_path": output_path,
                "delimiter": delimiter,
                "plot_snap_ogs_output": plot_snap_ogs_output,
                "plot_format": plot_format,
                "total_taxa": total_taxa,
            },
        )

    return {
        "status": "completed",
        "subgroup_counter": subgroup_counter,
        "subgroup_records": subgroup_records,
    }


def _coerce_float(value, fallback):
    if value is None or value == "":
        return fallback
    return float(value)


def _coerce_int(value, fallback):
    if value is None or value == "":
        return fallback
    return int(value)


def _execute_manifest_runs(config: dict):
    manifest_path = config["manifest"]
    output_root = config["output_path"]
    if not output_root.endswith("/"):
        output_root += "/"
    os.makedirs(output_root, exist_ok=True)

    delimiter = "\t" if manifest_path.endswith(".tsv") else ","

    summary_rows = []
    with open(manifest_path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        required = {"tree", "fasta"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise SystemExit("Manifest must include 'tree' and 'fasta' columns.")

        for idx, row in enumerate(reader, start=1):
            tree = row.get("tree")
            fasta = row.get("fasta")
            run_output_path = row.get("output_path") or output_root
            if not run_output_path.endswith("/"):
                run_output_path += "/"

            if row.get("id"):
                run_output_path = f"{run_output_path}{row['id']}/"

            os.makedirs(run_output_path, exist_ok=True)

            run_cfg = dict(config)
            run_cfg.update(
                {
                    "manifest": None,
                    "tree": tree,
                    "fasta": fasta,
                    "support": _coerce_float(row.get("support"), config["support"]),
                    "occupancy": _coerce_float(row.get("occupancy"), config["occupancy"]),
                    "occupancy_count": _coerce_int(row.get("occupancy_count"), config.get("occupancy_count")),
                    "occupancy_fraction": _coerce_float(row.get("occupancy_fraction"), config.get("occupancy_fraction")),
                    "rooted": _parse_bool(row.get("rooted"), config["rooted"]),
                    "snap_trees": _parse_bool(row.get("snap_trees"), config["snap_trees"]),
                    "report_inparalog_handling": _parse_bool(
                        row.get("report_inparalog_handling"), config["report_inparalog_handling"]
                    ),
                    "output_path": run_output_path,
                    "delimiter": row.get("delimiter") or config["delimiter"],
                }
            )

            inparalog_value = row.get("inparalog_to_keep")
            if inparalog_value:
                run_cfg["inparalog_to_keep"] = InparalogToKeep(inparalog_value)

            if run_cfg.get("occupancy_fraction") is not None:
                run_cfg["occupancy_mode"] = "fraction"
                unique_taxa = len(
                    {
                        rec.id.split(run_cfg["delimiter"], 1)[0]
                        for rec in SeqIO.parse(fasta, "fasta")
                    }
                )
                run_cfg["occupancy"] = max(
                    1,
                    int((run_cfg["occupancy_fraction"] * unique_taxa) + 0.999999),
                )
            elif run_cfg.get("occupancy_count") is not None:
                run_cfg["occupancy_mode"] = "count"
                run_cfg["occupancy"] = run_cfg["occupancy_count"]
            elif run_cfg.get("occupancy") is None:
                run_cfg["occupancy"] = determine_occupancy_threshold(
                    fasta, run_cfg["delimiter"]
                )

            run_cfg.pop("manifest", None)
            result = execute(**run_cfg)
            summary_rows.append(
                {
                    "row": idx,
                    "tree": tree,
                    "fasta": fasta,
                    "status": result.get("status", "completed"),
                    "subgroup_count": result.get("subgroup_counter", 0),
                    "output_path": run_output_path,
                }
            )

    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%SZ")
    summary_tsv = f"{output_root}manifest_summary_{stamp}.tsv"
    summary_json = f"{output_root}manifest_summary_{stamp}.json"

    with open(summary_tsv, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["row", "tree", "fasta", "status", "subgroup_count", "output_path"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    with open(summary_json, "w") as handle:
        json.dump(summary_rows, handle, indent=2)

    print(f"Manifest execution summary TSV: {summary_tsv}")
    print(f"Manifest execution summary JSON: {summary_json}")


def main(argv=None):
    """
    Function that parses and collects arguments
    """
    parser = create_parser()
    if argv is None:
        argv = sys.argv[1:]
    if len(argv) == 0:
        parser.print_help(sys.stderr)
        sys.exit(2)

    args = parser.parse_args(argv)
    config = process_args(args)

    if config.get("manifest"):
        _execute_manifest_runs(config)
    else:
        execute_config = dict(config)
        execute_config.pop("manifest", None)
        execute(**execute_config)


if __name__ == "__main__":
    main(sys.argv[1:])
