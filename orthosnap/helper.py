from collections import Counter
import copy
from enum import Enum
import re
import statistics as stat

from Bio import Phylo
from Bio import SeqIO
from Bio.Phylo.BaseTree import TreeMixin

class InparalogToKeep(Enum):
    shortest_seq_len = "shortest_seq_len"
    median_seq_len = "median_seq_len"
    longest_seq_len = "longest_seq_len"
    shortest_branch_len = "shortest_branch_len"
    median_branch_len = "median_branch_len"
    longest_branch_len = "longest_branch_len"

def collapse_low_support_bipartitions(newtree, support: float):
    """
    collapse bipartitions with support less than threshold
    """

    newtree.collapse_all(lambda c: c.confidence is not None and c.confidence < support)

    return newtree


def determine_if_dups_are_sister(subtree_tips: list):
    """
    determine if dups are sister to one another
    """
    # get first set of subtree tips
    first_set_of_subtree_tips = subtree_tips[0]
    # set if duplicate sequences are sister as True
    are_sisters = True
    # check if duplicate sequences are sister
    for set_of_subtree_tips in subtree_tips[1:]:
        if first_set_of_subtree_tips != set_of_subtree_tips:
            are_sisters = False
        if not are_sisters:
            break

    return are_sisters


def get_all_tips_and_taxa_names(tree):
    """
    get all taxa and tip names in a phylogeny

    return lists with information from each
    """
    taxa = []
    all_tips = []

    # loop through the tree and collect terminal names
    for term in tree.get_terminals():
        taxa_name = term.name[: term.name.index("|")]
        if taxa_name not in taxa:
            taxa.append(taxa_name)
        all_tips.append(term.name)

    return taxa, all_tips

def check_if_single_copy(taxa: list, all_tips: list):
    """
    check if the input phylogeny is already a single-copy tree
    """

    if len(taxa) == len(all_tips):
        print("Input phylogeny is already a single-copy orthogroup\nExiting now...")
        return True
    else:
        return False


def get_tips_and_taxa_names_and_taxa_counts_from_subtrees(inter):
    """
    get taxa, counts of each taxa, and all terminal names
    """
    taxa_from_terms = []
    terms = []
    # get taxa and terminal names from subtree
    for term in inter.get_terminals():
        taxa_from_terms.append(term.name.split("|", 1)[0])
        terms.append(term.name)
    # count number of taxa in subtree
    counts_of_taxa_from_terms = Counter(taxa_from_terms)
    counts = []
    # count number of times each taxon is present
    for count in counts_of_taxa_from_terms.values():
        counts.append(count)

    return taxa_from_terms, terms, counts_of_taxa_from_terms, counts


def get_subtree_tips(terms: list, name: str, tree):
    """
    get lists of subsubtrees from subtree
    """
    # get the duplicate sequences
    dups = [e for e in terms if e.startswith(name)]
    subtree_tips = []
    # for individual sequence among duplicate sequences
    for dup in dups:
        # create a copy of the tree
        temptree = copy.deepcopy(tree)
        # get the node path for the duplicate sequence
        node_path = temptree.get_path(dup)
        # for the terminals of the parent of the duplicate sequence
        # get the terminal names and append them to temp
        temp = []
        for term in node_path[-2].get_terminals():
            temp.append(term.name)
        subtree_tips.append(temp)

    # print(name)
    # import sys
    # sys.exit()

    return subtree_tips, dups


def handle_multi_copy_subtree(
    all_tips: list,
    terms: list,
    newtree,
    subgroup_counter: int,
    fasta: str,
    support: float,
    fasta_dict: dict,
    assigned_tips: list,
    counts_of_taxa_from_terms,
    tree,
    snap_trees: bool,
    inparalog_to_keep: InparalogToKeep,
):
    """
    handling case where subtree contains all single copy genes
    """
    # prune subtree to get subtree of interest
    newtree = prune_subtree(all_tips, terms, newtree)

    # collapse bipartition with low support
    newtree = collapse_low_support_bipartitions(newtree, support)

    # remove duplicate sequences if they are sister to one another
    # following the approach in PhyloTreePruner
    for name in counts_of_taxa_from_terms:
        # if the taxon is represented by more than one sequence
        if counts_of_taxa_from_terms[name] > 1:
            # get subtree tips
            subtree_tips, dups = get_subtree_tips(terms, name, tree)

            # check if subtrees are sister to one another
            are_sisters = determine_if_dups_are_sister(subtree_tips)

            # if duplicate sequences are sister, get the longest sequence
            if are_sisters:
                # trim short sequences and keep long sequences in newtree
                newtree, terms = inparalog_to_keep_determination(newtree, fasta_dict, dups, terms, inparalog_to_keep)

    # if the resulting subtree has only single copy genes
    # create a fasta file with sequences from tip labels
    _, _, _, counts = get_tips_and_taxa_names_and_taxa_counts_from_subtrees(newtree)
    if set(counts) == set([1]):
        (
            subgroup_counter,
            assigned_tips,
        ) = write_output_fasta_and_account_for_assigned_tips_single_copy_case(
            fasta,
            subgroup_counter,
            terms,
            fasta_dict,
            assigned_tips,
            snap_trees,
            newtree,
        )

    return subgroup_counter, assigned_tips


def handle_single_copy_subtree(
    all_tips: list,
    terms: list,
    newtree,
    subgroup_counter: int,
    fasta: str,
    support: float,
    fasta_dict: dict,
    assigned_tips: list,
    snap_trees: bool,
):
    """
    handling case where subtree contains all single copy genes
    """
    # prune subtree to get subtree of interest
    newtree = prune_subtree(all_tips, terms, newtree)

    # collapse bipartition with low support
    newtree = collapse_low_support_bipartitions(newtree, support)

    # add list of terms to assigned_tips list
    # and create subgroup fasta files
    (
        subgroup_counter,
        assigned_tips,
    ) = write_output_fasta_and_account_for_assigned_tips_single_copy_case(
        fasta, subgroup_counter, terms, fasta_dict, assigned_tips, snap_trees, newtree
    )

    return subgroup_counter, assigned_tips


def inparalog_to_keep_determination(
    newtree,
    fasta_dict: dict,
    dups: list,
    terms: list,
    inparalog_to_keep: InparalogToKeep,
):
    """
    remove_short_sequences_among_duplicates_that_are_sister
    """

    lengths = dict()
    # keep inparalog based on sequence length
    if inparalog_to_keep.value in ["shortest_seq_len", "median_seq_len", "longest_seq_len"]:
        for dup in dups:
            lengths[dup] = len(re.sub("-", "", str(fasta_dict[dup].seq)))
        # determine which sequences to keep
        if inparalog_to_keep.value == "shortest_seq_len":
            seq_to_keep = min(lengths, key=lengths.get)
        elif len(lengths) > 2 and inparalog_to_keep.value == "median_seq_len":
            median_len = stat.median(lengths, key=lengths.get)
            seq_to_keep = [key for key,value in lengths if value==median_len]
        elif len(lengths) == 2 and inparalog_to_keep.value == "median_seq_len":
            seq_to_keep = max(lengths, key=lengths.get)
        elif inparalog_to_keep.value == "longest_seq_len":
            seq_to_keep = max(lengths, key=lengths.get)
            
    # keep inparalog based on tip to root length
    else:
        for dup in dups:
            lengths[dup] = TreeMixin.distance(newtree, dup)
        if inparalog_to_keep.value == "shortest_branch_len":
            seq_to_keep = min(lengths, key=lengths.get)
        elif len(lengths) > 2 and inparalog_to_keep.value == "median_branch_len":
            median_len = stat.median(lengths, key=lengths.get)
            seq_to_keep = [key for key,value in lengths if value==median_len]
        elif len(lengths) == 2 and inparalog_to_keep.value == "median_branch_len":
            seq_to_keep = max(lengths, key=lengths.get)
        elif inparalog_to_keep.value == "longest_branch_len":
            seq_to_keep = max(lengths, key=lengths.get)

    # trim unwanted species-specific
    # paralogous sequences from the tree
    for seq_name, _ in lengths.items():
        if seq_name != seq_to_keep:
            newtree.prune(seq_name)
            terms.remove(seq_name)

    return newtree, terms


def prune_subtree(all_tips: list, terms: list, newtree):
    """
    prune tips not of interest for subtree
    """

    tips_to_prune = [i for i in all_tips + terms if i not in all_tips or i not in terms]

    for tip in tips_to_prune:
        newtree.prune(tip)

    return newtree


def read_input_files(tree: str, fasta: str, rooted: bool):
    """
    read input files and midpoint root tree
    """

    tree = Phylo.read(tree, "newick")

    if not rooted:
        tree.root_at_midpoint()

    fasta = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    return tree, fasta


def write_output_fasta_and_account_for_assigned_tips_single_copy_case(
    fasta: str,
    subgroup_counter: int,
    terms: list,
    fasta_dict: dict,
    assigned_tips: list,
    snap_tree: bool,
    newtree,
):

    # write output
    output_file_name = f"{fasta}.orthosnap.{subgroup_counter}.fa"
    with open(output_file_name, "w") as output_handle:
        for term in terms:
            SeqIO.write(fasta_dict[term], output_handle, "fasta")
            assigned_tips.append(term)

    if snap_tree:
        output_file_name = f"{fasta}.orthosnap.{subgroup_counter}.tre"
        Phylo.write(newtree, output_file_name, "newick")

    subgroup_counter += 1

    return subgroup_counter, assigned_tips
