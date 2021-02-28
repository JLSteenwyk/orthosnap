import logging
import os.path
import sys

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
    occupancy = args.occupancy if args.occupancy is not None else 0.5

    return dict(
        tree=tree,
        fasta=fasta,
        support=support,
        occupancy=occupancy
    )