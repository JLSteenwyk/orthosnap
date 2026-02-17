from collections import Counter
from enum import Enum
import re
import sys

from Bio import Phylo
from Bio import SeqIO
from Bio.Phylo.BaseTree import TreeMixin, Tree


class InparalogToKeep(Enum):
    shortest_seq_len = "shortest_seq_len"
    median_seq_len = "median_seq_len"
    longest_seq_len = "longest_seq_len"
    shortest_branch_len = "shortest_branch_len"
    median_branch_len = "median_branch_len"
    longest_branch_len = "longest_branch_len"


def clone_subtree_as_tree(subtree):
    """Clone a clade subtree into a standalone Bio.Phylo Tree."""

    def _clone_clade(clade):
        clone = clade.__class__()
        for key, value in clade.__dict__.items():
            if key != "clades":
                setattr(clone, key, value)
        clone.clades = []
        return clone

    root_clone = _clone_clade(subtree)
    stack = [(subtree, root_clone)]

    while stack:
        source, target = stack.pop()
        for child in source.clades:
            child_clone = _clone_clade(child)
            target.clades.append(child_clone)
            stack.append((child, child_clone))

    return Tree(root=root_clone)


def collapse_low_support_bipartitions(newtree, support: float):
    """
    collapse bipartitions with support less than threshold
    """

    newtree.collapse_all(lambda c: c.confidence is not None and c.confidence < support)

    return newtree


def determine_if_dups_are_sister(
    subtree_tips: list,
    clade_terminal_sets: set,
):
    """
    determine if dups are sister to one another
    """
    return frozenset(subtree_tips) in clade_terminal_sets


def build_clade_terminal_set_index(tree):
    """
    Build a set of terminal-name frozensets for every clade in a tree.
    """

    clade_to_terminals = dict()
    terminal_sets = set()

    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            names = frozenset([clade.name])
        else:
            merged = set()
            for child in clade.clades:
                merged.update(clade_to_terminals[child])
            names = frozenset(merged)

        clade_to_terminals[clade] = names
        terminal_sets.add(names)

    return terminal_sets


def update_clade_terminal_set_index_for_pruned_tips(
    clade_terminal_sets: set,
    pruned_tips: list,
):
    """
    Incrementally update terminal-set index after pruning one or more tips.
    """

    updated_sets = clade_terminal_sets
    for tip in pruned_tips:
        next_sets = set()
        for terms in updated_sets:
            if tip in terms:
                trimmed_terms = terms.difference([tip])
                if trimmed_terms:
                    next_sets.add(frozenset(trimmed_terms))
            else:
                next_sets.add(terms)
        updated_sets = next_sets

    return updated_sets


def get_all_tips_and_taxa_names(tree, delimiter: str):
    """
    get all taxa and tip names in a phylogeny

    return lists with information from each
    """
    taxa = []
    all_tips = []

    # loop through the tree and collect terminal names
    for term in tree.get_terminals():
        try:
            taxa_name = term.name[: term.name.index(delimiter)]
        except ValueError:
            print("\nERROR: Delimiter does not exist in FASTA headers.\nSpecify the delimiter using the -d argument.")
            sys.exit()
        if taxa_name not in taxa:
            taxa.append(taxa_name)
        all_tips.append(term.name)

    return taxa, all_tips


def build_tip_parent_lookup(tree):
    """Return a mapping from terminal name to its parent clade."""

    tip_parent = dict()

    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            if child.is_terminal() and child.name is not None:
                tip_parent[child.name] = parent

    return tip_parent


def check_if_single_copy(taxa: list, all_tips: list):
    """
    check if the input phylogeny is already a single-copy tree
    """

    if len(taxa) == len(all_tips):
        print("Input phylogeny is already a single-copy orthogroup\nExiting now...")
        return True
    else:
        return False


def get_tips_and_taxa_names_and_taxa_counts_from_subtrees(inter, delimiter: str):
    """
    get taxa, counts of each taxa, and all terminal names
    """
    taxa_from_terms = []
    terms = []
    # get taxa and terminal names from subtree
    for term in inter.get_terminals():
        taxa_from_terms.append(term.name.split(delimiter, 1)[0])
        terms.append(term.name)
    # count number of taxa in subtree
    counts_of_taxa_from_terms = Counter(taxa_from_terms)
    # count number of times each taxon is present
    counts = list(counts_of_taxa_from_terms.values())

    return taxa_from_terms, terms, counts_of_taxa_from_terms, counts


def build_subtree_taxa_cache(tree, delimiter: str):
    """
    Cache subtree term and taxa-count data for each internal clade.
    """

    term_cache = dict()
    taxa_count_cache = dict()
    subtree_cache = dict()

    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            term_name = clade.name
            term_cache[clade] = [term_name]
            taxa_count_cache[clade] = Counter([term_name.split(delimiter, 1)[0]])
            continue

        terms = []
        counts_of_taxa_from_terms = Counter()

        for child in clade.clades:
            terms.extend(term_cache[child])
            counts_of_taxa_from_terms.update(taxa_count_cache[child])

        term_cache[clade] = terms
        taxa_count_cache[clade] = counts_of_taxa_from_terms
        subtree_cache[clade] = (
            terms,
            set(terms),
            counts_of_taxa_from_terms,
            list(counts_of_taxa_from_terms.values()),
        )

    return subtree_cache


def get_subtree_tips(terms: list, name: str):
    """
    get lists of subsubtrees from subtree
    """
    # get the duplicate sequences
    dups = [e for e in terms if e.startswith(name)]
    # This helper is only used for duplicate discovery.
    # Keep return shape for compatibility with older callers/tests.
    return [], dups


def handle_multi_copy_subtree(
    subtree,
    terms: list,
    subgroup_counter: int,
    fasta: str,
    support: float,
    fasta_dict: dict,
    assigned_tips: set,
    counts_of_taxa_from_terms,
    snap_trees: bool,
    inparalog_to_keep: InparalogToKeep,
    output_path: str,
    inparalog_handling: dict,
    inparalog_handling_summary: dict,
    report_inparalog_handling: bool,
    delimiter: str,
):
    """
    handling case where subtree contains all single copy genes
    """
    newtree = clone_subtree_as_tree(subtree)

    # collapse bipartition with low support
    newtree = collapse_low_support_bipartitions(newtree, support)
    clade_terminal_sets = build_clade_terminal_set_index(newtree)

    # remove duplicate sequences if they are sister to one another
    # following the approach in PhyloTreePruner
    for name in counts_of_taxa_from_terms:
        # if the taxon is represented by more than one sequence
        if counts_of_taxa_from_terms[name] > 1:
            # get subtree tips
            _, dups = get_subtree_tips(terms, name)
            if not dups:
                continue

            # check if subtrees are sister to one another
            are_sisters = determine_if_dups_are_sister(
                dups, clade_terminal_sets
            )

            # if duplicate sequences are sister, get the longest sequence
            if are_sisters:
                # trim short sequences and keep long sequences in newtree
                newtree, terms, inparalog_handling, pruned_tips = \
                    inparalog_to_keep_determination(
                        newtree, fasta_dict, dups, terms,
                        inparalog_to_keep, inparalog_handling
                    )
                clade_terminal_sets = update_clade_terminal_set_index_for_pruned_tips(
                    clade_terminal_sets, pruned_tips
                )

    # if the resulting subtree has only single copy genes
    # create a fasta file with sequences from tip labels
    _, _, _, counts = \
        get_tips_and_taxa_names_and_taxa_counts_from_subtrees(newtree, delimiter)

    if set(counts) == set([1]):
        (
            subgroup_counter,
            assigned_tips,
            inparalog_handling_summary
        ) = write_output_fasta_and_account_for_assigned_tips_single_copy_case(
            fasta,
            subgroup_counter,
            terms,
            fasta_dict,
            assigned_tips,
            snap_trees,
            newtree,
            output_path,
            inparalog_handling,
            inparalog_handling_summary,
            report_inparalog_handling,
        )

    return \
        subgroup_counter, assigned_tips, \
        inparalog_handling, inparalog_handling_summary


def handle_single_copy_subtree(
    subtree,
    terms: list,
    subgroup_counter: int,
    fasta: str,
    support: float,
    fasta_dict: dict,
    assigned_tips: set,
    snap_trees: bool,
    output_path: str,
    inparalog_handling: dict,
    inparalog_handling_summary: dict,
    report_inparalog_handling: bool,
):
    """
    handling case where subtree contains all single copy genes
    """
    newtree = clone_subtree_as_tree(subtree)

    # collapse bipartition with low support
    newtree = collapse_low_support_bipartitions(newtree, support)

    # add list of terms to assigned_tips list
    # and create subgroup fasta files
    (
        subgroup_counter,
        assigned_tips,
        inparalog_handling_summary
    ) = write_output_fasta_and_account_for_assigned_tips_single_copy_case(
        fasta,
        subgroup_counter,
        terms,
        fasta_dict,
        assigned_tips,
        snap_trees,
        newtree,
        output_path,
        inparalog_handling,
        inparalog_handling_summary,
        report_inparalog_handling,
    )

    return \
        subgroup_counter, assigned_tips, \
        inparalog_handling, inparalog_handling_summary


def inparalog_to_keep_determination(
    newtree,
    fasta_dict: dict,
    dups: list,
    terms: list,
    inparalog_to_keep: InparalogToKeep,
    inparalog_handling: dict,
):
    """
    remove_short_sequences_among_duplicates_that_are_sister
    """
    lengths = dict()
    pruned_tips = []
    def _select_median_key(values: dict):
        """Select deterministic median key by value, tie-broken by key."""
        sorted_items = sorted(values.items(), key=lambda item: (item[1], item[0]))
        return sorted_items[len(sorted_items) // 2][0]

    # keep inparalog based on sequence length
    if inparalog_to_keep.value in [
        "shortest_seq_len",
        "median_seq_len",
        "longest_seq_len",
    ]:
        for dup in dups:
            lengths[dup] = len(re.sub("-", "", str(fasta_dict[dup].seq)))
        # determine which sequences to keep
        if inparalog_to_keep.value == "shortest_seq_len":
            seq_to_keep = min(lengths, key=lengths.get)
        elif len(lengths) > 2 and inparalog_to_keep.value == "median_seq_len":
            seq_to_keep = _select_median_key(lengths)
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
        elif len(lengths) > 2 and \
                inparalog_to_keep.value == "median_branch_len":
            seq_to_keep = _select_median_key(lengths)
        elif len(lengths) == 2 and \
                inparalog_to_keep.value == "median_branch_len":
            seq_to_keep = max(lengths, key=lengths.get)
        elif inparalog_to_keep.value == "longest_branch_len":
            seq_to_keep = max(lengths, key=lengths.get)

    # trim unwanted species-specific
    # paralogous sequences from the tree
    for seq_name, _ in lengths.items():
        if seq_name != seq_to_keep:
            newtree.prune(seq_name)
            terms.remove(seq_name)
            pruned_tips.append(seq_name)

    inparalog_handling[seq_to_keep] = [dup for dup in dups if dup != seq_to_keep]

    return newtree, terms, inparalog_handling, pruned_tips


def prune_subtree(all_tips: list, terms: list, newtree):
    """
    prune tips not of interest for subtree
    """

    tips_to_prune = [
        i for i in all_tips + terms if i not in all_tips or i not in terms
    ]

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
    assigned_tips: set,
    snap_tree: bool,
    newtree,
    output_path: str,
    inparalog_handling: dict,
    inparalog_handling_summary: dict,
    report_inparalog_handling: bool,
):
    # write output
    fasta_path_stripped = re.sub("^.*/", "", fasta)
    output_file_name = (
        f"{output_path}/{fasta_path_stripped}.orthosnap.{subgroup_counter}.fa"
    )
    with open(output_file_name, "w") as output_handle:
        for term in terms:
            SeqIO.write(fasta_dict[term], output_handle, "fasta")
            assigned_tips.add(term)

    if snap_tree:
        output_file_name = (
            f"{output_path}/{fasta_path_stripped}.orthosnap.{subgroup_counter}.tre"
        )
        Phylo.write(newtree, output_file_name, "newick")

    if report_inparalog_handling:
        write_summary_file_with_inparalog_handling(
            inparalog_handling,
            fasta,
            output_path,
            subgroup_counter,
            terms,
        )
    subgroup_counter += 1

    return subgroup_counter, assigned_tips, inparalog_handling_summary


def write_summary_file_with_inparalog_handling(
        inparalog_handling: dict,
        fasta: str,
        output_path: str,
        subgroup_count: int,
        kept_terms: list = None,
):
    res_arr = []

    in_file_handle = re.sub("^.*/", "", fasta)

    for k, v in inparalog_handling.items():
        temp = []
        temp.append(in_file_handle+".orthosnap."+str(subgroup_count))
        temp.append(k)
        temp.append(';'.join(v))
        res_arr.append(temp)
    inparalog_report_output_name = in_file_handle + ".inparalog_report.txt"

    if kept_terms is None:
        kept_terms = []

    kept_terms_set = set(kept_terms)
    with open(f"{output_path}{inparalog_report_output_name}", "a") as file:
        for i in res_arr:
            if i[1] in kept_terms_set:
                file.writelines('\t'.join(i) + '\n')
    pruned_tips = []
