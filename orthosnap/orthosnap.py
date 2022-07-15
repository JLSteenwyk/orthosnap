#!/usr/bin/env python

import copy
import sys
import time

from tqdm import tqdm

from .args_processing import process_args
from .helper import (
    check_if_single_copy,
    get_all_tips_and_taxa_names,
    get_tips_and_taxa_names_and_taxa_counts_from_subtrees,
    handle_single_copy_subtree,
    handle_multi_copy_subtree,
    read_input_files,
)
from .helper import InparalogToKeep
from .parser import create_parser
from .writer import write_user_args, write_output_stats


def execute(
    tree: str,
    fasta: str,
    support: float,
    occupancy: float,
    rooted: bool,
    snap_trees: bool,
    inparalog_to_keep: InparalogToKeep,
):
    """
    Master execute Function
    -----------------------
    This function executes the main functions and calls other subfunctions
    """

    # write user args to stdout
    write_user_args(
        tree, fasta,
        support, occupancy,
        rooted, snap_trees,
        inparalog_to_keep
    )

    # create start time logger
    start_time = time.time()

    # read input files and midpoint root tree
    tree, fasta_dict = read_input_files(tree, fasta, rooted)

    # get list of all tip names and taxa names
    taxa, all_tips = get_all_tips_and_taxa_names(tree)

    # check if the inputted phylogeny is already a single-copy tree.
    # if it is, exit
    if check_if_single_copy(taxa, all_tips):
        sys.exit()

    # loop through tree, but skip the root (hence [1:])
    # keep tabs of terms that have already been assigned
    # to a subgroup as well as a counter for that subgroup
    assigned_tips = []
    subgroup_counter = 0
    
    for inter in tqdm(tree.get_nonterminals()[1:]):
        (
            _,
            terms,
            counts_of_taxa_from_terms,
            counts,
        ) = get_tips_and_taxa_names_and_taxa_counts_from_subtrees(inter)

        # create a copy of the input tree
        newtree = copy.deepcopy(tree)

        # if a sufficient number of taxa are represented, examine the subtree
        if len(counts_of_taxa_from_terms) >= occupancy:
            # if each taxon is represented by one sequence and
            # the tips have not been assigned to a suborthogroup
            # prune tips not part of the subtree of interest
            if (
                set([1]) == set(counts)
                and len(list(set(terms) & set(assigned_tips))) == 0
            ):
                subgroup_counter, assigned_tips = handle_single_copy_subtree(
                    all_tips,
                    terms,
                    newtree,
                    subgroup_counter,
                    fasta,
                    support,
                    fasta_dict,
                    assigned_tips,
                    snap_trees,
                )

            # if any taxon is represented by more than one sequence and
            # the tips have not been assigned to a suborthogroup
            # prune tips not part of the subtree of interest
            elif len(list(set(terms) & set(assigned_tips))) == 0:
                subgroup_counter, assigned_tips = handle_multi_copy_subtree(
                    all_tips,
                    terms,
                    newtree,
                    subgroup_counter,
                    fasta,
                    support,
                    fasta_dict,
                    assigned_tips,
                    counts_of_taxa_from_terms,
                    tree,
                    snap_trees,
                    inparalog_to_keep,
                )

    write_output_stats(fasta, subgroup_counter, start_time, snap_trees)


def main(argv=None):
    """
    Function that parses and collects arguments
    """
    # parse and assign arguments
    parser = create_parser()
    args = parser.parse_args()

    # pass to master execute function
    execute(**process_args(args))


if __name__ == "__main__":
    main(sys.argv[1:])
