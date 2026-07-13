from collections import Counter
from collections.abc import Mapping
from enum import Enum
from bisect import bisect_left, bisect_right, insort
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


def _clone_clade_without_children(clade):
    clone = clade.__class__()
    for key, value in clade.__dict__.items():
        if key != "clades":
            setattr(clone, key, value)
    clone.clades = []
    return clone


def clone_subtree_as_tree(subtree):
    """Clone a clade subtree into a standalone Bio.Phylo Tree."""
    root_clone = _clone_clade_without_children(subtree)
    stack = [(subtree, root_clone)]

    while stack:
        source, target = stack.pop()
        for child in source.clades:
            child_clone = _clone_clade_without_children(child)
            target.clades.append(child_clone)
            stack.append((child, child_clone))

    return Tree(root=root_clone)


def clone_induced_tree(tree, keep_tips: set):
    """Clone the minimal tree induced by ``keep_tips`` without a full copy."""
    cloned_clades = {}
    stack = [(tree.root, False)]

    while stack:
        clade, visited = stack.pop()
        if clade.is_terminal():
            cloned_clades[clade] = (
                _clone_clade_without_children(clade)
                if clade.name in keep_tips
                else None
            )
            continue
        if not visited:
            stack.append((clade, True))
            stack.extend((child, False) for child in reversed(clade.clades))
            continue

        selected_children = []
        for child in clade.clades:
            child_clone = cloned_clades.pop(child)
            if child_clone is not None:
                selected_children.append(child_clone)
        if not selected_children:
            cloned_clades[clade] = None
        elif len(selected_children) == 1 and len(clade.clades) > 1:
            child = selected_children[0]
            if child.branch_length is not None:
                child.branch_length += clade.branch_length or 0.0
            cloned_clades[clade] = child
        else:
            clone = _clone_clade_without_children(clade)
            clone.clades = selected_children
            cloned_clades[clade] = clone

    root_clone = cloned_clades[tree.root]
    if root_clone is None:
        raise ValueError("none of the requested tips exist in the reference tree")

    tree_clone = tree.__class__()
    for key, value in tree.__dict__.items():
        if key != "root":
            setattr(tree_clone, key, value)
    tree_clone.root = root_clone
    return tree_clone


def collapse_low_support_bipartitions(newtree, support: float):
    """
    collapse bipartitions with support less than threshold
    """

    newtree.collapse_all(lambda c: c.confidence is not None and c.confidence < support)

    return newtree


def determine_if_dups_are_sister(
    subtree_tips: list,
    clade_terminal_index: dict,
):
    """
    determine if dups are sister to one another
    """
    tip_bits = clade_terminal_index["tip_bits"]
    duplicate_mask = 0
    for tip in subtree_tips:
        bit = tip_bits.get(tip)
        if bit is None:
            return False
        duplicate_mask |= bit
    return duplicate_mask in clade_terminal_index["masks"]


def build_clade_terminal_set_index(tree):
    """
    Build a bitmask index for terminal membership across all clades.
    """

    terminal_names = sorted(
        term.name for term in tree.get_terminals() if term.name is not None
    )
    tip_bits = {name: (1 << idx) for idx, name in enumerate(terminal_names)}
    clade_to_mask = dict()
    masks = set()

    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            mask = tip_bits[clade.name]
        else:
            merged_mask = 0
            for child in clade.clades:
                merged_mask |= clade_to_mask[child]
            mask = merged_mask

        clade_to_mask[clade] = mask
        masks.add(mask)
    return {"masks": masks, "tip_bits": tip_bits}


def update_clade_terminal_set_index_for_pruned_tips(
    clade_terminal_index: dict,
    pruned_tips: list,
):
    """
    Incrementally update bitmask index after pruning one or more tips.
    """

    if not pruned_tips:
        return clade_terminal_index

    prune_mask = 0
    tip_bits = clade_terminal_index["tip_bits"]
    for tip in pruned_tips:
        bit = tip_bits.get(tip)
        if bit is not None:
            prune_mask |= bit

    if prune_mask == 0:
        return clade_terminal_index

    updated_masks = set()
    for mask in clade_terminal_index["masks"]:
        trimmed_mask = mask & ~prune_mask
        if trimmed_mask:
            updated_masks.add(trimmed_mask)

    return {"masks": updated_masks, "tip_bits": tip_bits}


def build_terminal_parent_maps(tree):
    """
    Build terminal-node and parent lookups for fast terminal pruning.
    """

    parent_lookup = dict()
    terminal_lookup = dict()

    stack = [tree.root]
    while stack:
        parent = stack.pop()
        for child in parent.clades:
            parent_lookup[child] = parent
            if child.is_terminal() and child.name is not None:
                terminal_lookup[child.name] = child
            else:
                stack.append(child)

    return terminal_lookup, parent_lookup


def prune_terminal_fast(tree, terminal_name: str, terminal_lookup: dict, parent_lookup: dict):
    """
    Prune one terminal by name using parent links to avoid full-tree traversals.
    """

    target = terminal_lookup.pop(terminal_name, None)
    if target is None:
        raise ValueError("can't find a matching target below this root")

    parent = parent_lookup.get(target)
    if parent is None:
        raise ValueError("can't find parent for matching target")

    parent.clades.remove(target)
    parent_lookup.pop(target, None)

    if len(parent.clades) == 1:
        child = parent.clades[0]
        if child.branch_length is not None:
            child.branch_length += parent.branch_length or 0.0

        grandparent = parent_lookup.get(parent)
        if grandparent is None:
            tree.root = child
            parent_lookup.pop(parent, None)
            if not child.is_terminal():
                for grandchild in child.clades:
                    parent_lookup[grandchild] = child
        else:
            idx = grandparent.clades.index(parent)
            grandparent.clades[idx] = child
            parent_lookup[child] = grandparent
            parent_lookup.pop(parent, None)

    return tree


def get_all_tips_and_taxa_names(tree, delimiter: str):
    """
    get all taxa and tip names in a phylogeny

    return lists with information from each
    """
    taxa = []
    seen_taxa = set()
    all_tips = []

    # loop through the tree and collect terminal names
    for term in tree.get_terminals():
        try:
            taxa_name = term.name[: term.name.index(delimiter)]
        except ValueError:
            print("\nERROR: Delimiter does not exist in FASTA headers.\nSpecify the delimiter using the -d argument.")
            sys.exit()
        if taxa_name not in seen_taxa:
            seen_taxa.add(taxa_name)
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


class SubtreeTaxaCache(Mapping):
    """Compact, lazily materialized terminal and taxon data for clades."""

    def __init__(self, tree, delimiter: str):
        self.delimiter = delimiter
        self.terminal_names = []
        self.taxon_names = []
        self.terminal_indices = {}
        self.terminal_clades = {}
        self.internal_clades = []
        self.parent_lookup = {}
        self._intervals = {}
        self._summaries = {}

        stack = [tree.root]
        while stack:
            clade = stack.pop()
            if clade.is_terminal():
                index = len(self.terminal_names)
                self.terminal_names.append(clade.name)
                self.taxon_names.append(clade.name.split(delimiter, 1)[0])
                self.terminal_indices[clade.name] = index
                self.terminal_clades[clade.name] = clade
                self._intervals[clade] = (index, index + 1)
                continue

            self.internal_clades.append(clade)
            for child in clade.clades:
                self.parent_lookup[child] = clade
            stack.extend(reversed(clade.clades))

        taxa_sets = {}
        stack = [(tree.root, False)]
        while stack:
            clade, visited = stack.pop()
            if clade.is_terminal():
                index = self._intervals[clade][0]
                taxa_sets[clade] = {self.taxon_names[index]}
                continue
            if not visited:
                stack.append((clade, True))
                stack.extend((child, False) for child in reversed(clade.clades))
                continue

            child_sets = [taxa_sets.pop(child) for child in clade.clades]
            largest_index = max(range(len(child_sets)), key=lambda i: len(child_sets[i]))
            taxa = child_sets[largest_index]
            for index, child_taxa in enumerate(child_sets):
                if index != largest_index:
                    taxa.update(child_taxa)

            start = self._intervals[clade.clades[0]][0]
            end = self._intervals[clade.clades[-1]][1]
            terminal_count = end - start
            self._intervals[clade] = (start, end)
            self._summaries[clade] = (start, end, len(taxa), terminal_count)
            taxa_sets[clade] = taxa

        taxa_sets.clear()

    def __getitem__(self, clade):
        terms, terms_set, counts_of_taxa = self.details(clade)
        return terms, terms_set, counts_of_taxa, list(counts_of_taxa.values())

    def __iter__(self):
        return iter(self._summaries)

    def __len__(self):
        return len(self._summaries)

    def summary(self, clade):
        """Return start, end, unique-taxon count, and terminal count."""
        return self._summaries[clade]

    def details(self, clade):
        """Materialize the legacy cache values for one clade only."""
        start, end, _, _ = self._summaries[clade]
        terms = self.terminal_names[start:end]
        counts_of_taxa = Counter(self.taxon_names[start:end])
        return terms, set(terms), counts_of_taxa

    def taxon_tip_groups(self, clade):
        """Return terminal names grouped by taxon in terminal traversal order."""
        start, end, _, _ = self._summaries[clade]
        groups = {}
        for index in range(start, end):
            groups.setdefault(self.taxon_names[index], []).append(
                self.terminal_names[index]
            )
        return groups

    def duplicates_are_sister(
        self,
        candidate,
        duplicate_tips: list,
        pruned_indices: list,
        support: float,
    ):
        """Test sisterhood after support collapse and earlier logical pruning."""
        positions = [self.terminal_indices[name] for name in duplicate_tips]
        first_position = min(positions)
        last_position = max(positions)
        clade = self.terminal_clades[duplicate_tips[0]]

        while True:
            start, end = self._intervals[clade]
            contains_duplicates = start <= first_position and last_position < end
            survives_collapse = (
                clade is candidate
                or clade.confidence is None
                or clade.confidence >= support
            )
            if contains_duplicates and survives_collapse:
                break
            if clade is candidate:
                return False
            clade = self.parent_lookup.get(clade)
            if clade is None:
                return False

        pruned_in_clade = bisect_right(pruned_indices, end - 1) - bisect_left(
            pruned_indices, start
        )
        return (end - start - pruned_in_clade) == len(duplicate_tips)


def build_subtree_taxa_cache(tree, delimiter: str):
    """
    Cache subtree term and taxa-count data for each internal clade.
    """

    return SubtreeTaxaCache(tree, delimiter)


def get_subtree_tips(terms: list, name: str, delimiter: str):
    """
    get lists of subsubtrees from subtree
    """
    # get the duplicate sequences
    dups = [e for e in terms if e.split(delimiter, 1)[0] == name]
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
    subgroup_records: list = None,
    write_outputs: bool = True,
    subtree_cache: SubtreeTaxaCache = None,
    taxon_tip_groups: dict = None,
    sequence_lengths: dict = None,
):
    """
    handling case where subtree contains all single copy genes
    """
    if subtree_cache is not None:
        return _handle_multi_copy_subtree_indexed(
            subtree=subtree,
            terms=terms,
            subgroup_counter=subgroup_counter,
            fasta=fasta,
            support=support,
            fasta_dict=fasta_dict,
            assigned_tips=assigned_tips,
            counts_of_taxa_from_terms=counts_of_taxa_from_terms,
            snap_trees=snap_trees,
            inparalog_to_keep=inparalog_to_keep,
            output_path=output_path,
            inparalog_handling=inparalog_handling,
            inparalog_handling_summary=inparalog_handling_summary,
            report_inparalog_handling=report_inparalog_handling,
            subgroup_records=subgroup_records,
            write_outputs=write_outputs,
            subtree_cache=subtree_cache,
            taxon_tip_groups=taxon_tip_groups,
            sequence_lengths=sequence_lengths,
        )

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
            _, dups = get_subtree_tips(terms, name, delimiter)
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
            subgroup_records,
            write_outputs,
        )

    return \
        subgroup_counter, assigned_tips, \
        inparalog_handling, inparalog_handling_summary


def _select_inparalog_to_keep(
    tree,
    fasta_dict: dict,
    dups: list,
    inparalog_to_keep: InparalogToKeep,
    sequence_lengths: dict = None,
):
    """Select one inparalog with the historical ordering and tie behavior."""
    values = {}

    if inparalog_to_keep.value in [
        "shortest_seq_len",
        "median_seq_len",
        "longest_seq_len",
    ]:
        if sequence_lengths is None:
            sequence_lengths = {}
        for dup in dups:
            if dup not in sequence_lengths:
                sequence_lengths[dup] = len(re.sub("-", "", str(fasta_dict[dup].seq)))
            values[dup] = sequence_lengths[dup]
    else:
        for dup in dups:
            values[dup] = TreeMixin.distance(tree, dup)

    if inparalog_to_keep.value in ["shortest_seq_len", "shortest_branch_len"]:
        return min(values, key=values.get)
    if inparalog_to_keep.value in ["longest_seq_len", "longest_branch_len"]:
        return max(values, key=values.get)
    if len(values) == 2:
        return max(values, key=values.get)

    sorted_items = sorted(values.items(), key=lambda item: (item[1], item[0]))
    return sorted_items[len(sorted_items) // 2][0]


def _handle_multi_copy_subtree_indexed(
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
    subgroup_records: list,
    write_outputs: bool,
    subtree_cache: SubtreeTaxaCache,
    taxon_tip_groups: dict,
    sequence_lengths: dict,
):
    """Process a multi-copy candidate using the shared topology index."""
    if taxon_tip_groups is None:
        taxon_tip_groups = subtree_cache.taxon_tip_groups(subtree)

    pruned_indices = []
    pruned_tips = []
    resolved_duplicate_taxa = set()
    distance_tree = Tree(root=subtree)

    for name, count in counts_of_taxa_from_terms.items():
        if count <= 1:
            continue
        dups = taxon_tip_groups[name]
        if not subtree_cache.duplicates_are_sister(
            subtree,
            dups,
            pruned_indices,
            support,
        ):
            continue

        seq_to_keep = _select_inparalog_to_keep(
            distance_tree,
            fasta_dict,
            dups,
            inparalog_to_keep,
            sequence_lengths,
        )
        trimmed = [dup for dup in dups if dup != seq_to_keep]
        inparalog_handling[seq_to_keep] = trimmed
        resolved_duplicate_taxa.add(name)
        pruned_tips.extend(trimmed)
        for tip in trimmed:
            insort(pruned_indices, subtree_cache.terminal_indices[tip])

    duplicate_taxa = {
        name for name, count in counts_of_taxa_from_terms.items() if count > 1
    }
    if resolved_duplicate_taxa == duplicate_taxa:
        pruned_tip_set = set(pruned_tips)
        kept_terms = [term for term in terms if term not in pruned_tip_set]
        newtree = None

        if write_outputs and snap_trees:
            newtree = clone_subtree_as_tree(subtree)
            collapse_low_support_bipartitions(newtree, support)
            terminal_lookup, parent_lookup = build_terminal_parent_maps(newtree)
            for tip in pruned_tips:
                prune_terminal_fast(newtree, tip, terminal_lookup, parent_lookup)

        (
            subgroup_counter,
            assigned_tips,
            inparalog_handling_summary,
        ) = write_output_fasta_and_account_for_assigned_tips_single_copy_case(
            fasta,
            subgroup_counter,
            kept_terms,
            fasta_dict,
            assigned_tips,
            snap_trees,
            newtree,
            output_path,
            inparalog_handling,
            inparalog_handling_summary,
            report_inparalog_handling,
            subgroup_records,
            write_outputs,
        )

    return (
        subgroup_counter,
        assigned_tips,
        inparalog_handling,
        inparalog_handling_summary,
    )


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
    subgroup_records: list = None,
    write_outputs: bool = True,
):
    """
    handling case where subtree contains all single copy genes
    """
    newtree = None
    if write_outputs and snap_trees:
        newtree = clone_subtree_as_tree(subtree)
        collapse_low_support_bipartitions(newtree, support)

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
        subgroup_records,
        write_outputs,
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
    pruned_tips = []
    seq_to_keep = _select_inparalog_to_keep(
        newtree,
        fasta_dict,
        dups,
        inparalog_to_keep,
    )

    # trim unwanted species-specific
    # paralogous sequences from the tree
    for seq_name in dups:
        if seq_name != seq_to_keep:
            pruned_tips.append(seq_name)

    if pruned_tips:
        terminal_lookup, parent_lookup = build_terminal_parent_maps(newtree)
        for tip_name in pruned_tips:
            prune_terminal_fast(newtree, tip_name, terminal_lookup, parent_lookup)
        pruned_tip_set = set(pruned_tips)
        terms = [term for term in terms if term not in pruned_tip_set]

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
    subgroup_records: list = None,
    write_outputs: bool = True,
):
    fasta_path_stripped = re.sub("^.*/", "", fasta)

    for term in terms:
        assigned_tips.add(term)

    if write_outputs:
        output_file_name = (
            f"{output_path}/{fasta_path_stripped}.orthosnap.{subgroup_counter}.fa"
        )
        with open(output_file_name, "w") as output_handle:
            SeqIO.write((fasta_dict[term] for term in terms), output_handle, "fasta")

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

    if subgroup_records is not None:
        subgroup_records.append(
            {"subgroup_id": subgroup_counter, "tips": list(terms)}
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
