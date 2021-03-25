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

    if occupancy <= 0:
        logger.warning("Occupancy threshold must be greater than 0.")
        sys.exit()

    if support > 100 or support < 0:
        logger.warning("Support threshold must range from 0 to 100.")
        logger.warning("I recommend using a threshold of 80 for UFBoot.")
        logger.warning("and 70 for classic bootstrap support.")
        sys.exit()

    return dict(tree=tree, fasta=fasta, support=support, occupancy=occupancy)
