import logging
import math
import os.path
import sys

from Bio import SeqIO

logger = logging.getLogger(__name__)


def process_args(args) -> dict:
    """
    Process args from argparser and set defaults
    """
    tree = args.tree
    fasta = args.fasta

    if not os.path.isfile(tree):
        logger.warning("Input tree does not exist")
        sys.exit()
    if not os.path.isfile(fasta):
        logger.warning("Input fasta does not exist")
        sys.exit()

    # assign optional arguments
    support = args.support if args.support is not None else 80
    occupancy = (
        args.occupancy
        if args.occupancy is not None
        else determine_occupancy_threshold(fasta)
    )

    if occupancy <= 0:
        logger.warning("Occupancy threshold must be greater than 0.")
        sys.exit()

    if support > 100 or support < 0:
        logger.warning("Support threshold must range from 0 to 100.")
        logger.warning("I recommend using a threshold of 80 for ultrafast")
        logger.warning("bootstrap approximations (also known as UFBoot)")
        logger.warning("and 70 for classic bootstrap support.")
        sys.exit()

    return dict(tree=tree, fasta=fasta, support=support, occupancy=occupancy)


def determine_occupancy_threshold(fasta: str) -> int:
    fasta = SeqIO.parse(fasta, "fasta")
    unique_names = []

    for i in fasta:
        i = i.id.split("|", 1)[0]
        if i not in unique_names:
            unique_names.append(i)

    occupancy_threshold = proper_round(len(unique_names) / 2)

    return occupancy_threshold


def proper_round(num):
    if num - math.floor(num) < 0.5:
        return(math.floor(num))
    else:
        return(math.ceil(num))
