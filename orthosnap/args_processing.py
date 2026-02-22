import logging
import math
import os.path
import re
import sys
from Bio import SeqIO

from .helper import InparalogToKeep

logger = logging.getLogger(__name__)


def process_args(args) -> dict:
    """
    Process args from argparser and set defaults
    """
    tree = args.tree
    fasta = args.fasta
    manifest = getattr(args, "manifest", None)

    delimiter = args.delimiter if args.delimiter is not None else "|"

    if manifest:
        if not os.path.isfile(manifest):
            logger.warning("Manifest file does not exist")
            sys.exit()
    else:
        if tree is None or fasta is None:
            logger.warning("Provide both --tree and --fasta, or provide --manifest.")
            sys.exit()
        if not os.path.isfile(tree):
            logger.warning("Input tree does not exist")
            sys.exit()
        if not os.path.isfile(fasta):
            logger.warning("Input fasta does not exist")
            sys.exit()

    # assign optional arguments
    support = args.support if args.support is not None else 80
    occupancy_count = getattr(args, "occupancy_count", None)
    occupancy_fraction = getattr(args, "occupancy_fraction", None)
    occupancy = args.occupancy

    occupancy_modes = [occupancy is not None, occupancy_count is not None, occupancy_fraction is not None]
    if sum(occupancy_modes) > 1:
        logger.warning("Use only one occupancy mode: -o/--occupancy, --occupancy-count, or --occupancy-fraction.")
        sys.exit()

    if occupancy_fraction is not None and (occupancy_fraction <= 0 or occupancy_fraction > 1):
        logger.warning("Occupancy fraction must be in the range (0, 1].")
        sys.exit()

    if occupancy_count is not None and occupancy_count <= 0:
        logger.warning("Occupancy count must be greater than 0.")
        sys.exit()

    if support > 100 or support < 0:
        logger.warning("Support threshold must range from 0 to 100.")
        logger.warning("I recommend using a threshold of 80 for ultrafast")
        logger.warning("bootstrap approximations (also known as UFBoot)")
        logger.warning("and 70 for classic bootstrap support.")
        sys.exit()

    rooted = args.rooted
    snap_trees = args.snap_trees
    report_inparalog_handling = args.report_inparalog_handling
    plot_snap_ogs_output = getattr(args, "plot_snap_ogs", False)
    raw_plot_format = getattr(args, "plot_format", None)
    plot_format = raw_plot_format if raw_plot_format is not None else "png"
    validate_only = getattr(args, "validate_only", False)
    resume = getattr(args, "resume", False)
    structured_output = getattr(args, "structured_output", False)
    bootstrap_trees = getattr(args, "bootstrap_trees", None)
    consensus_min_frequency = getattr(args, "consensus_min_frequency", None)
    consensus_trees = getattr(args, "consensus_trees", False)

    if consensus_min_frequency is None:
        consensus_min_frequency = 0.5
    if consensus_min_frequency <= 0 or consensus_min_frequency > 1:
        logger.warning("Consensus minimum frequency must be in the range (0, 1].")
        sys.exit()

    if bootstrap_trees is not None and not os.path.isfile(bootstrap_trees):
        logger.warning("Bootstrap tree list file does not exist.")
        sys.exit()

    if args.output_path:
        output_path = args.output_path
        if not output_path.endswith("/"):
            output_path = output_path + "/"
    else:
        if fasta is None:
            output_path = "./"
        else:
            output_path = re.sub("/[^/]+$", "", fasta)
            if output_path == fasta:
                output_path = "./"
            elif not output_path.endswith("/"):
                output_path = output_path + "/"

    if args.inparalog_to_keep:
        inparalog_to_keep = InparalogToKeep(args.inparalog_to_keep)
    else:
        inparalog_to_keep = InparalogToKeep.longest_seq_len

    occupancy_mode = "legacy"
    if occupancy_count is not None:
        occupancy_mode = "count"
    elif occupancy_fraction is not None:
        occupancy_mode = "fraction"

    if fasta is not None:
        total_taxa = count_unique_taxa_in_fasta(fasta, delimiter)
    else:
        total_taxa = None

    resolved_occupancy = None
    if occupancy_mode == "count":
        resolved_occupancy = int(occupancy_count)
    elif occupancy_mode == "fraction":
        if total_taxa is None:
            logger.warning("Cannot resolve occupancy fraction without a FASTA input.")
            sys.exit()
        resolved_occupancy = max(1, int(math.ceil(occupancy_fraction * total_taxa)))
    elif occupancy is not None:
        resolved_occupancy = occupancy
    elif fasta is not None:
        resolved_occupancy = determine_occupancy_threshold(fasta, delimiter)

    if resolved_occupancy is not None and resolved_occupancy <= 0:
        logger.warning("Occupancy threshold must be greater than 0.")
        sys.exit()

    return dict(
        tree=tree,
        fasta=fasta,
        support=support,
        occupancy=resolved_occupancy,
        occupancy_mode=occupancy_mode,
        occupancy_count=occupancy_count,
        occupancy_fraction=occupancy_fraction,
        rooted=rooted,
        snap_trees=snap_trees,
        inparalog_to_keep=inparalog_to_keep,
        report_inparalog_handling=report_inparalog_handling,
        output_path=output_path,
        delimiter=delimiter,
        plot_snap_ogs_output=plot_snap_ogs_output,
        plot_format=plot_format,
        manifest=manifest,
        validate_only=validate_only,
        resume=resume,
        structured_output=structured_output,
        bootstrap_trees=bootstrap_trees,
        consensus_min_frequency=consensus_min_frequency,
        consensus_trees=consensus_trees,
        total_taxa=total_taxa,
    )


def determine_occupancy_threshold(fasta: str, delimiter: str) -> int:
    fasta = SeqIO.parse(fasta, "fasta")
    unique_names = []

    for i in fasta:
        i = i.id.split(delimiter, 1)[0]
        if i not in unique_names:
            unique_names.append(i)

    occupancy_threshold = proper_round(len(unique_names) / 2)

    return occupancy_threshold


def count_unique_taxa_in_fasta(fasta: str, delimiter: str) -> int:
    """Count unique taxon labels in a FASTA using the configured delimiter."""
    fasta_handle = SeqIO.parse(fasta, "fasta")
    unique_names = set()

    for seq in fasta_handle:
        seq_id = seq.id
        if delimiter in seq_id:
            unique_names.add(seq_id.split(delimiter, 1)[0])
        else:
            unique_names.add(seq_id)
    return len(unique_names)


def proper_round(num):
    if num - math.floor(num) < 0.5:
        return math.floor(num)
    else:
        return math.ceil(num)
