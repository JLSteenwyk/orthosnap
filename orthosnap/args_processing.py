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
    occupancy = args.occupancy if args.occupancy is not None else 1.0

    if occupancy > 1 or occupancy < 0:
        logger.warning("Occupancy threshold must range from 0 to 1.")
        logger.warning("I recommend using a high threshold close to one.")
        logger.warning("For example, 0.5, 0.75, or 1.")
        sys.exit()

    return dict(tree=tree, fasta=fasta, support=support, occupancy=occupancy)
