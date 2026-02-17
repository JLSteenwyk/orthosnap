from collections import Counter
import copy
from io import StringIO
from pathlib import Path
import pytest

from Bio import Phylo
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from orthosnap.helper import (
    build_subtree_taxa_cache,
    build_clade_terminal_set_index,
    check_if_single_copy,
    clone_subtree_as_tree,
    collapse_low_support_bipartitions,
    determine_if_dups_are_sister,
    get_all_tips_and_taxa_names,
    get_tips_and_taxa_names_and_taxa_counts_from_subtrees,
    get_subtree_tips,
    inparalog_to_keep_determination,
    build_terminal_parent_maps,
    prune_terminal_fast,
    prune_subtree,
    read_input_files,
    update_clade_terminal_set_index_for_pruned_tips,
    write_output_fasta_and_account_for_assigned_tips_single_copy_case,
)
from orthosnap.helper import InparalogToKeep

here = Path(__file__)


def _index_to_terminal_sets(index):
    bit_to_tip = {bit: tip for tip, bit in index["tip_bits"].items()}
    terminal_sets = set()

    for mask in index["masks"]:
        members = set()
        remaining = mask
        while remaining:
            lowest_bit = remaining & -remaining
            members.add(bit_to_tip[lowest_bit])
            remaining ^= lowest_bit
        terminal_sets.add(frozenset(members))

    return terminal_sets


class TestCollapseLowSupportBipartitions(object):
    def test_collapses_bipartitions(self):
        ## setup
        tree = Phylo.read(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            "newick",
        )
        tree_copy = copy.deepcopy(tree)
        tree.collapse_all(lambda c: c.confidence is not None and c.confidence < 80)

        ## execution
        support = 80
        res = collapse_low_support_bipartitions(tree_copy, support)

        ## check results
        for term0, term1 in zip(tree.get_terminals(), res.get_terminals()):
            assert term0.name == term1.name
            assert term0.branch_length == term1.branch_length


class TestCloneSubtreeAsTree(object):
    def test_clone_subtree_as_tree_independent_copy(self):
        tree = Phylo.read(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            "newick",
        )
        subtree = tree.common_ancestor(
            [
                "Aspergillus_fumigatus_Af293|EAL84942.1",
                "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_06972-RA",
            ]
        )

        cloned_tree = clone_subtree_as_tree(subtree)
        clone_terms = {term.name for term in cloned_tree.get_terminals()}
        source_terms = {term.name for term in subtree.get_terminals()}

        assert clone_terms == source_terms

        original_children = len(subtree.clades)
        cloned_tree.root.clades.pop()

        assert len(subtree.clades) == original_children
        assert len(cloned_tree.root.clades) == original_children - 1


class TestDetermineIfDupsAreSister(object):
    def test_determine_if_dups_are_sister_true(self):
        tree = Phylo.read(StringIO("((sp|a,sp|b),sp|c);"), "newick")
        clade_terminal_sets = build_clade_terminal_set_index(tree)

        assert determine_if_dups_are_sister(["sp|a", "sp|b"], clade_terminal_sets)

    def test_determine_if_dups_are_sister_false(self):
        tree = Phylo.read(StringIO("((sp|a,sp|b),sp|c);"), "newick")
        clade_terminal_sets = build_clade_terminal_set_index(tree)

        assert not determine_if_dups_are_sister(
            ["sp|a", "sp|c"], clade_terminal_sets
        )


class TestUpdateCladeTerminalSetIndexForPrunedTips(object):
    def test_incremental_index_matches_rebuild(self):
        tree = Phylo.read(StringIO("(((sp|a,sp|b),sp|c),sp|d);"), "newick")
        before = build_clade_terminal_set_index(tree)

        tree.prune("sp|d")
        tree.prune("sp|c")
        rebuilt = build_clade_terminal_set_index(tree)

        updated = update_clade_terminal_set_index_for_pruned_tips(
            before, ["sp|d", "sp|c"]
        )

        assert _index_to_terminal_sets(updated) == _index_to_terminal_sets(rebuilt)

    def test_repeated_pruned_tip_is_idempotent(self):
        tree = Phylo.read(StringIO("(((sp|a,sp|b),sp|c),sp|d);"), "newick")
        before = build_clade_terminal_set_index(tree)

        once = update_clade_terminal_set_index_for_pruned_tips(
            copy.deepcopy(before), ["sp|d"]
        )
        repeated = update_clade_terminal_set_index_for_pruned_tips(
            copy.deepcopy(before), ["sp|d", "sp|d"]
        )

        assert repeated == once

    def test_missing_pruned_tip_does_not_change_index(self):
        tree = Phylo.read(StringIO("(((sp|a,sp|b),sp|c),sp|d);"), "newick")
        before = build_clade_terminal_set_index(tree)

        updated = update_clade_terminal_set_index_for_pruned_tips(
            copy.deepcopy(before), ["sp|not_present"]
        )

        assert updated == before

    def test_prune_to_singleton_matches_rebuild(self):
        tree = Phylo.read(StringIO("((sp|a,sp|b),sp|c);"), "newick")
        before = build_clade_terminal_set_index(tree)

        tree.prune("sp|b")
        tree.prune("sp|c")
        rebuilt = build_clade_terminal_set_index(tree)

        updated = update_clade_terminal_set_index_for_pruned_tips(
            before, ["sp|b", "sp|c"]
        )

        assert _index_to_terminal_sets(updated) == _index_to_terminal_sets(rebuilt)


# class TestDetermineIfDupsAreSister(object):
#     def test_determine_if_dups_are_sister_true(self):
#         ## setup
#         subtree_tips = [
#             [
#                 "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#                 "Aspergillus_fumigatus_Af293|EAL85095.2",
#             ],
#             [
#                 "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#                 "Aspergillus_fumigatus_Af293|EAL85095.2",
#             ],
#         ]
#         tree = Phylo.read(
#             f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
#             "newick",
#         )
#         tree = tree.common_ancestor(subtree_tips[0])
#         expected_res = True

#         ## execution
#         res = determine_if_dups_are_sister(subtree_tips, tree)

#         ## check results
#         assert res == expected_res

#     def test_determine_if_dups_are_sister_false(self):
#         ## setup
#         subtree_tips = [
#             [
#                 "Aspergillus_niger_CBS_513.88|XP_001391581.1",
#                 "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
#                 "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
#                 "Aspergillus_awamori_IFM_58123|GCB17486.1",
#             ],
#             [
#                 "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#                 "Aspergillus_awamori_IFM_58123|GCB19008.1",
#             ],
#         ]
#         expected_res = False

#         ## execution
#         res = determine_if_dups_are_sister(subtree_tips)

#         ## check results
#         assert res == expected_res


class TestGetAllTipsAndTaxaNames(object):
    def test_get_all_tips_and_taxa_names0(self):
        ## setup
        tree = Phylo.read(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            "newick",
        )
        tree.root_at_midpoint()
        expected_all_tips = [
            "Aspergillus_fumigatus_Af293|EAL84942.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_10038-RA",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_06972-RA",
            "Aspergillus_fumigatus_Af293|EAL85274.1",
            "Aspergillus_fischeri_NRRL181|XP_001262055.1",
            "Aspergillus_niger_CBS_513.88|XP_001398067.1",
            "Aspergillus_awamori_IFM_58123|GCB24888.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_09150-RA",
            "Aspergillus_fumigatus_Af293|EAL84262.1",
            "Aspergillus_fischeri_NRRL181|XP_001267441.1",
            "Aspergillus_niger_CBS_513.88|XP_001400898.1",
            "Aspergillus_awamori_IFM_58123|GCB19337.1",
            "Aspergillus_niger_CBS_513.88|XP_001397349.2",
            "Aspergillus_awamori_IFM_58123|GCB25420.1",
            "Aspergillus_fumigatus_Af293|EAL87116.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_03570-RA",
            "Aspergillus_fischeri_NRRL181|XP_001265570.1",
            "Aspergillus_niger_CBS_513.88|XP_001390417.1",
            "Aspergillus_niger_CBS_513.88|XP_001390415.2",
            "Aspergillus_awamori_IFM_58123|GCB24160.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_08472-RA",
            "Aspergillus_fumigatus_Af293|EAL94045.1",
            "Aspergillus_fischeri_NRRL181|XP_001261225.1",
            "Aspergillus_niger_CBS_513.88|XP_001402298.1",
            "Aspergillus_awamori_IFM_58123|GCB23392.1",
            "Aspergillus_fumigatus_Af293|EAL93843.2",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_04905-RA",
            "Aspergillus_fischeri_NRRL181|XP_001261009.1",
            "Aspergillus_niger_CBS_513.88|XP_001397083.2",
            "Aspergillus_awamori_IFM_58123|GCB27653.1",
            "Aspergillus_niger_CBS_513.88|XP_001391581.1",
            "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
            "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
            "Aspergillus_awamori_IFM_58123|GCB17486.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
            "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
            "Aspergillus_fumigatus_Af293|EAL85095.2",
            "Aspergillus_fischeri_NRRL181|XP_001261692.1",
            "Aspergillus_niger_CBS_513.88|XP_001401336.1",
            "Aspergillus_awamori_IFM_58123|GCB19008.1",
            "Aspergillus_awamori_IFM_58123|GCB22318.1",
        ]
        expected_taxa = [
            "Aspergillus_fumigatus_Af293",
            "Aspergillus_oerlinghausenensis_CBS139183",
            "Aspergillus_fischeri_NRRL181",
            "Aspergillus_niger_CBS_513.88",
            "Aspergillus_awamori_IFM_58123",
        ]

        ## execution
        taxa, all_tips = get_all_tips_and_taxa_names(tree, "|")

        ## check results
        assert set(taxa) == set(expected_taxa)
        assert set(all_tips) == set(expected_all_tips)


class TestGetTipsAndTaxaNamesAndTaxaCountsFromSubtrees(object):
    def test_get_all_tips_and_taxa_names1(self):
        ## setup
        tree = Phylo.read(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            "newick",
        )
        tree.root_at_midpoint()
        expected_taxa_from_terms = [
            "Aspergillus_fumigatus_Af293",
            "Aspergillus_oerlinghausenensis_CBS139183",
            "Aspergillus_oerlinghausenensis_CBS139183",
            "Aspergillus_fumigatus_Af293",
            "Aspergillus_fischeri_NRRL181",
            "Aspergillus_niger_CBS_513.88",
            "Aspergillus_awamori_IFM_58123",
            "Aspergillus_oerlinghausenensis_CBS139183",
            "Aspergillus_fumigatus_Af293",
            "Aspergillus_fischeri_NRRL181",
            "Aspergillus_niger_CBS_513.88",
            "Aspergillus_awamori_IFM_58123",
            "Aspergillus_niger_CBS_513.88",
            "Aspergillus_awamori_IFM_58123",
            "Aspergillus_fumigatus_Af293",
            "Aspergillus_oerlinghausenensis_CBS139183",
            "Aspergillus_fischeri_NRRL181",
            "Aspergillus_niger_CBS_513.88",
            "Aspergillus_niger_CBS_513.88",
            "Aspergillus_awamori_IFM_58123",
            "Aspergillus_oerlinghausenensis_CBS139183",
            "Aspergillus_fumigatus_Af293",
            "Aspergillus_fischeri_NRRL181",
            "Aspergillus_niger_CBS_513.88",
            "Aspergillus_awamori_IFM_58123",
            "Aspergillus_fumigatus_Af293",
            "Aspergillus_oerlinghausenensis_CBS139183",
            "Aspergillus_fischeri_NRRL181",
            "Aspergillus_niger_CBS_513.88",
            "Aspergillus_awamori_IFM_58123",
            "Aspergillus_niger_CBS_513.88",
            "Aspergillus_awamori_IFM_58123",
            "Aspergillus_awamori_IFM_58123",
            "Aspergillus_awamori_IFM_58123",
            "Aspergillus_oerlinghausenensis_CBS139183",
            "Aspergillus_oerlinghausenensis_CBS139183",
            "Aspergillus_fumigatus_Af293",
            "Aspergillus_fumigatus_Af293",
            "Aspergillus_fischeri_NRRL181",
            "Aspergillus_niger_CBS_513.88",
            "Aspergillus_awamori_IFM_58123",
        ]
        expected_terms = [
            "Aspergillus_awamori_IFM_58123|GCB22318.1",
            "Aspergillus_fumigatus_Af293|EAL93843.2",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_04905-RA",
            "Aspergillus_fischeri_NRRL181|XP_001261009.1",
            "Aspergillus_niger_CBS_513.88|XP_001397083.2",
            "Aspergillus_awamori_IFM_58123|GCB27653.1",
            "Aspergillus_niger_CBS_513.88|XP_001397349.2",
            "Aspergillus_awamori_IFM_58123|GCB25420.1",
            "Aspergillus_fumigatus_Af293|EAL87116.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_03570-RA",
            "Aspergillus_fischeri_NRRL181|XP_001265570.1",
            "Aspergillus_niger_CBS_513.88|XP_001390417.1",
            "Aspergillus_niger_CBS_513.88|XP_001390415.2",
            "Aspergillus_awamori_IFM_58123|GCB24160.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_08472-RA",
            "Aspergillus_fumigatus_Af293|EAL94045.1",
            "Aspergillus_fischeri_NRRL181|XP_001261225.1",
            "Aspergillus_niger_CBS_513.88|XP_001402298.1",
            "Aspergillus_awamori_IFM_58123|GCB23392.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_10038-RA",
            "Aspergillus_fumigatus_Af293|EAL84942.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_06972-RA",
            "Aspergillus_fumigatus_Af293|EAL85274.1",
            "Aspergillus_fischeri_NRRL181|XP_001262055.1",
            "Aspergillus_niger_CBS_513.88|XP_001398067.1",
            "Aspergillus_awamori_IFM_58123|GCB24888.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_09150-RA",
            "Aspergillus_fumigatus_Af293|EAL84262.1",
            "Aspergillus_fischeri_NRRL181|XP_001267441.1",
            "Aspergillus_niger_CBS_513.88|XP_001400898.1",
            "Aspergillus_awamori_IFM_58123|GCB19337.1",
        ]
        expected_counts_of_taxa_from_terms = Counter(
            {
                "Aspergillus_awamori_IFM_58123": 7,
                "Aspergillus_niger_CBS_513.88": 7,
                "Aspergillus_fumigatus_Af293": 6,
                "Aspergillus_oerlinghausenensis_CBS139183": 6,
                "Aspergillus_fischeri_NRRL181": 5,
            }
        )
        expected_counts = [7, 6, 6, 5, 7]

        ## execution
        for inter in tree.get_nonterminals()[1:]:
            (
                taxa_from_terms,
                terms,
                counts_of_taxa_from_terms,
                counts,
            ) = get_tips_and_taxa_names_and_taxa_counts_from_subtrees(inter, "|")
            break

        ## check results
        assert set(taxa_from_terms) == set(expected_taxa_from_terms)
        assert set(terms) == set(expected_terms)
        assert counts_of_taxa_from_terms == expected_counts_of_taxa_from_terms
        assert counts == expected_counts


class TestBuildSubtreeTaxaCache(object):
    def test_cache_matches_subtree_function(self):
        tree = Phylo.read(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            "newick",
        )
        tree.root_at_midpoint()
        cache = build_subtree_taxa_cache(tree, "|")

        for inter in tree.get_nonterminals()[1:]:
            (
                _,
                terms,
                counts_of_taxa_from_terms,
                counts,
            ) = get_tips_and_taxa_names_and_taxa_counts_from_subtrees(inter, "|")
            (
                cached_terms,
                cached_terms_set,
                cached_counts_of_taxa_from_terms,
                cached_counts,
            ) = cache[inter]

            assert cached_terms == terms
            assert cached_terms_set == set(terms)
            assert cached_counts_of_taxa_from_terms == counts_of_taxa_from_terms
            assert cached_counts == counts


class TestGetSubtreeTips(object):
    def test_taxon_prefix_collision_is_not_matched(self):
        terms = ["sp1|geneA", "sp10|geneB", "sp1|geneC"]
        _, dups = get_subtree_tips(terms, "sp1", "|")
        assert dups == ["sp1|geneA", "sp1|geneC"]


class TestFastTerminalPrune(object):
    def test_fast_prune_matches_sequential_prune_for_terminals(self):
        tree_fast = Phylo.read(StringIO("(((a:1,b:1):1,c:1):1,d:1);"), "newick")
        tree_seq = Phylo.read(StringIO("(((a:1,b:1):1,c:1):1,d:1);"), "newick")

        terminal_lookup, parent_lookup = build_terminal_parent_maps(tree_fast)
        prune_terminal_fast(tree_fast, "b", terminal_lookup, parent_lookup)
        prune_terminal_fast(tree_fast, "d", terminal_lookup, parent_lookup)
        tree_seq.prune("b")
        tree_seq.prune("d")

        assert tree_fast.format("newick") == tree_seq.format("newick")


# class TestGetSubtreeTips(object):
#     def test_get_subtree_tips(self):
#         ## setup
#         tree = Phylo.read(
#             f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
#             "newick",
#         )
#         tree.root_at_midpoint()
#         expected_subtree_tips = [
#             [
#                 "Aspergillus_fumigatus_Af293|EAL84942.1",
#                 "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_10038-RA",
#             ],
#             [
#                 "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_06972-RA",
#                 "Aspergillus_fumigatus_Af293|EAL85274.1",
#             ],
#             [
#                 "Aspergillus_fumigatus_Af293|EAL84262.1",
#                 "Aspergillus_fischeri_NRRL181|XP_001267441.1",
#             ],
#             [
#                 "Aspergillus_fumigatus_Af293|EAL87116.1",
#                 "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_03570-RA",
#                 "Aspergillus_fischeri_NRRL181|XP_001265570.1",
#             ],
#             [
#                 "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_08472-RA",
#                 "Aspergillus_fumigatus_Af293|EAL94045.1",
#             ],
#             [
#                 "Aspergillus_fumigatus_Af293|EAL93843.2",
#                 "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_04905-RA",
#                 "Aspergillus_fischeri_NRRL181|XP_001261009.1",
#             ],
#             [
#                 "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#                 "Aspergillus_fumigatus_Af293|EAL85095.2",
#             ],
#             [
#                 "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#                 "Aspergillus_fumigatus_Af293|EAL85095.2",
#             ],
#         ]
#         expected_dups = [
#             "Aspergillus_fumigatus_Af293|EAL84942.1",
#             "Aspergillus_fumigatus_Af293|EAL85274.1",
#             "Aspergillus_fumigatus_Af293|EAL84262.1",
#             "Aspergillus_fumigatus_Af293|EAL87116.1",
#             "Aspergillus_fumigatus_Af293|EAL94045.1",
#             "Aspergillus_fumigatus_Af293|EAL93843.2",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2",
#         ]

#         ## execution
#         terms = [
#             "Aspergillus_fumigatus_Af293|EAL84942.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_10038-RA",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_06972-RA",
#             "Aspergillus_fumigatus_Af293|EAL85274.1",
#             "Aspergillus_fischeri_NRRL181|XP_001262055.1",
#             "Aspergillus_niger_CBS_513.88|XP_001398067.1",
#             "Aspergillus_awamori_IFM_58123|GCB24888.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_09150-RA",
#             "Aspergillus_fumigatus_Af293|EAL84262.1",
#             "Aspergillus_fischeri_NRRL181|XP_001267441.1",
#             "Aspergillus_niger_CBS_513.88|XP_001400898.1",
#             "Aspergillus_awamori_IFM_58123|GCB19337.1",
#             "Aspergillus_niger_CBS_513.88|XP_001397349.2",
#             "Aspergillus_awamori_IFM_58123|GCB25420.1",
#             "Aspergillus_fumigatus_Af293|EAL87116.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_03570-RA",
#             "Aspergillus_fischeri_NRRL181|XP_001265570.1",
#             "Aspergillus_niger_CBS_513.88|XP_001390417.1",
#             "Aspergillus_niger_CBS_513.88|XP_001390415.2",
#             "Aspergillus_awamori_IFM_58123|GCB24160.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_08472-RA",
#             "Aspergillus_fumigatus_Af293|EAL94045.1",
#             "Aspergillus_fischeri_NRRL181|XP_001261225.1",
#             "Aspergillus_niger_CBS_513.88|XP_001402298.1",
#             "Aspergillus_awamori_IFM_58123|GCB23392.1",
#             "Aspergillus_fumigatus_Af293|EAL93843.2",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_04905-RA",
#             "Aspergillus_fischeri_NRRL181|XP_001261009.1",
#             "Aspergillus_niger_CBS_513.88|XP_001397083.2",
#             "Aspergillus_awamori_IFM_58123|GCB27653.1",
#             "Aspergillus_niger_CBS_513.88|XP_001391581.1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2",
#             "Aspergillus_fischeri_NRRL181|XP_001261692.1",
#             "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#             "Aspergillus_awamori_IFM_58123|GCB19008.1",
#         ]

#         for inter in tree.get_nonterminals()[1:]:
#             subtree_tips, dups = get_subtree_tips(
#                 terms, "Aspergillus_awamori_IFM_58123", inter
#             )
#             break

#         ## check results
#         assert subtree_tips == expected_subtree_tips
#         assert dups == expected_dups


# class TestInparalogToKeepDetermination(object):
#     def test_inparalog_to_keep_determination_longest_seq_len(self):
#         ## setup
#         tree = Phylo.read(
#             f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
#             "newick",
#         )
#         tree.root_at_midpoint()
#         fasta_dict = SeqIO.to_dict(
#             SeqIO.parse(
#                 f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
#                 "fasta",
#             )
#         )
#         expected_terms = [
#             "Aspergillus_niger_CBS_513.88|XP_001391581.1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2",
#             "Aspergillus_fischeri_NRRL181|XP_001261692.1",
#             "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#             "Aspergillus_awamori_IFM_58123|GCB19008.1",
#         ]

#         ## execution
#         dups = [
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#         ]
#         terms = [
#             "Aspergillus_niger_CBS_513.88|XP_001391581.1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2",
#             "Aspergillus_fischeri_NRRL181|XP_001261692.1",
#             "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#             "Aspergillus_awamori_IFM_58123|GCB19008.1",
#         ]

#         for inter in tree.get_nonterminals()[31:]:
#             (_, terms) = inparalog_to_keep_determination(
#                 inter,
#                 fasta_dict,
#                 dups,
#                 terms,
#                 InparalogToKeep.longest_seq_len,
#             )
#             break

#         ## check results
#         assert terms == expected_terms

#     def test_inparalog_to_keep_determination_shortest_seq_len(self):
#         ## setup
#         tree = Phylo.read(
#             f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
#             "newick",
#         )
#         tree.root_at_midpoint()
#         fasta_dict = SeqIO.to_dict(
#             SeqIO.parse(
#                 f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
#                 "fasta",
#             )
#         )
#         expected_terms = [
#             "Aspergillus_niger_CBS_513.88|XP_001391581.1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2",
#             "Aspergillus_fischeri_NRRL181|XP_001261692.1",
#             "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#             "Aspergillus_awamori_IFM_58123|GCB19008.1",
#         ]

#         ## execution
#         dups = [
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#         ]
#         terms = [
#             "Aspergillus_niger_CBS_513.88|XP_001391581.1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2",
#             "Aspergillus_fischeri_NRRL181|XP_001261692.1",
#             "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#             "Aspergillus_awamori_IFM_58123|GCB19008.1",
#         ]

#         for inter in tree.get_nonterminals()[31:]:
#             (_, terms) = inparalog_to_keep_determination(
#                 inter,
#                 fasta_dict,
#                 dups,
#                 terms,
#                 InparalogToKeep.shortest_seq_len,
#             )
#             break

#         ## check results
#         assert terms == expected_terms

#     def test_inparalog_to_keep_determination_median_seq_len(self):
#         ## setup
#         tree = Phylo.read(
#             f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
#             "newick",
#         )
#         tree.root_at_midpoint()
#         fasta_dict = SeqIO.to_dict(
#             SeqIO.parse(
#                 f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
#                 "fasta",
#             )
#         )
#         expected_terms = [
#             "Aspergillus_niger_CBS_513.88|XP_001391581.1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2",
#             "Aspergillus_fischeri_NRRL181|XP_001261692.1",
#             "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#             "Aspergillus_awamori_IFM_58123|GCB19008.1",
#         ]

#         ## execution
#         dups = [
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#         ]
#         terms = [
#             "Aspergillus_niger_CBS_513.88|XP_001391581.1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2",
#             "Aspergillus_fischeri_NRRL181|XP_001261692.1",
#             "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#             "Aspergillus_awamori_IFM_58123|GCB19008.1",
#         ]

#         for inter in tree.get_nonterminals()[31:]:
#             (_, terms) = inparalog_to_keep_determination(
#                 inter,
#                 fasta_dict,
#                 dups,
#                 terms,
#                 InparalogToKeep.median_seq_len,
#             )
#             break

#         ## check results
#         assert terms == expected_terms

#     def test_inparalog_to_keep_determination_longest_branch_len(self):
#         ## setup
#         tree = Phylo.read(
#             f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
#             "newick",
#         )
#         tree.root_at_midpoint()
#         fasta_dict = SeqIO.to_dict(
#             SeqIO.parse(
#                 f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
#                 "fasta",
#             )
#         )
#         expected_terms = [
#             "Aspergillus_niger_CBS_513.88|XP_001391581.1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2",
#             "Aspergillus_fischeri_NRRL181|XP_001261692.1",
#             "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#             "Aspergillus_awamori_IFM_58123|GCB19008.1",
#         ]

#         ## execution
#         dups = [
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#         ]
#         terms = [
#             "Aspergillus_niger_CBS_513.88|XP_001391581.1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2",
#             "Aspergillus_fischeri_NRRL181|XP_001261692.1",
#             "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#             "Aspergillus_awamori_IFM_58123|GCB19008.1",
#         ]

#         for inter in tree.get_nonterminals()[31:]:
#             (_, terms) = inparalog_to_keep_determination(
#                 inter,
#                 fasta_dict,
#                 dups,
#                 terms,
#                 InparalogToKeep.longest_branch_len,
#             )
#             break

#         ## check results
#         assert terms == expected_terms

#     def test_inparalog_to_keep_determination_shortest_branch_len(self):
#         ## setup
#         tree = Phylo.read(
#             f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
#             "newick",
#         )
#         tree.root_at_midpoint()
#         fasta_dict = SeqIO.to_dict(
#             SeqIO.parse(
#                 f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
#                 "fasta",
#             )
#         )
#         expected_terms = [
#             "Aspergillus_niger_CBS_513.88|XP_001391581.1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2",
#             "Aspergillus_fischeri_NRRL181|XP_001261692.1",
#             "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#             "Aspergillus_awamori_IFM_58123|GCB19008.1",
#         ]

#         ## execution
#         dups = [
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#         ]
#         terms = [
#             "Aspergillus_niger_CBS_513.88|XP_001391581.1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2",
#             "Aspergillus_fischeri_NRRL181|XP_001261692.1",
#             "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#             "Aspergillus_awamori_IFM_58123|GCB19008.1",
#         ]

#         for inter in tree.get_nonterminals()[31:]:
#             (_, terms) = inparalog_to_keep_determination(
#                 inter,
#                 fasta_dict,
#                 dups,
#                 terms,
#                 InparalogToKeep.shortest_branch_len,
#             )
#             break

#         ## check results
#         assert terms == expected_terms

#     def test_inparalog_to_keep_determination_median_branch_len(self):
#         ## setup
#         tree = Phylo.read(
#             f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
#             "newick",
#         )
#         tree.root_at_midpoint()
#         fasta_dict = SeqIO.to_dict(
#             SeqIO.parse(
#                 f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
#                 "fasta",
#             )
#         )
#         expected_terms = [
#             "Aspergillus_niger_CBS_513.88|XP_001391581.1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2",
#             "Aspergillus_fischeri_NRRL181|XP_001261692.1",
#             "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#             "Aspergillus_awamori_IFM_58123|GCB19008.1",
#         ]

#         ## execution
#         dups = [
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#         ]
#         terms = [
#             "Aspergillus_niger_CBS_513.88|XP_001391581.1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
#             "Aspergillus_awamori_IFM_58123|GCB17486.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2",
#             "Aspergillus_fischeri_NRRL181|XP_001261692.1",
#             "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#             "Aspergillus_awamori_IFM_58123|GCB19008.1",
#         ]

#         for inter in tree.get_nonterminals()[31:]:
#             (_, terms) = inparalog_to_keep_determination(
#                 inter,
#                 fasta_dict,
#                 dups,
#                 terms,
#                 InparalogToKeep.median_branch_len,
#             )
#             break

#         ## check results
#         assert terms == expected_terms


class TestPruneSubtree(object):
    def test_prune_subtree(self):
        ## setup
        tree = Phylo.read(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            "newick",
        )
        tree.root_at_midpoint()
        newtree = copy.deepcopy(tree)
        tips_to_prune = ["Aspergillus_awamori_IFM_58123|GCB22318.1"]
        for tip in tips_to_prune:
            tree.prune(tip)

        ## execution
        all_tips = [
            "Aspergillus_fumigatus_Af293|EAL84942.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_10038-RA",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_06972-RA",
            "Aspergillus_fumigatus_Af293|EAL85274.1",
            "Aspergillus_fischeri_NRRL181|XP_001262055.1",
            "Aspergillus_niger_CBS_513.88|XP_001398067.1",
            "Aspergillus_awamori_IFM_58123|GCB24888.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_09150-RA",
            "Aspergillus_fumigatus_Af293|EAL84262.1",
            "Aspergillus_fischeri_NRRL181|XP_001267441.1",
            "Aspergillus_niger_CBS_513.88|XP_001400898.1",
            "Aspergillus_awamori_IFM_58123|GCB19337.1",
            "Aspergillus_niger_CBS_513.88|XP_001397349.2",
            "Aspergillus_awamori_IFM_58123|GCB25420.1",
            "Aspergillus_fumigatus_Af293|EAL87116.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_03570-RA",
            "Aspergillus_fischeri_NRRL181|XP_001265570.1",
            "Aspergillus_niger_CBS_513.88|XP_001390417.1",
            "Aspergillus_niger_CBS_513.88|XP_001390415.2",
            "Aspergillus_awamori_IFM_58123|GCB24160.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_08472-RA",
            "Aspergillus_fumigatus_Af293|EAL94045.1",
            "Aspergillus_fischeri_NRRL181|XP_001261225.1",
            "Aspergillus_niger_CBS_513.88|XP_001402298.1",
            "Aspergillus_awamori_IFM_58123|GCB23392.1",
            "Aspergillus_fumigatus_Af293|EAL93843.2",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_04905-RA",
            "Aspergillus_fischeri_NRRL181|XP_001261009.1",
            "Aspergillus_niger_CBS_513.88|XP_001397083.2",
            "Aspergillus_awamori_IFM_58123|GCB27653.1",
            "Aspergillus_niger_CBS_513.88|XP_001391581.1",
            "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
            "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
            "Aspergillus_awamori_IFM_58123|GCB17486.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
            "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
            "Aspergillus_fumigatus_Af293|EAL85095.2",
            "Aspergillus_fischeri_NRRL181|XP_001261692.1",
            "Aspergillus_niger_CBS_513.88|XP_001401336.1",
            "Aspergillus_awamori_IFM_58123|GCB19008.1",
            "Aspergillus_awamori_IFM_58123|GCB22318.1",
        ]
        terms = [
            "Aspergillus_fumigatus_Af293|EAL84942.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_10038-RA",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_06972-RA",
            "Aspergillus_fumigatus_Af293|EAL85274.1",
            "Aspergillus_fischeri_NRRL181|XP_001262055.1",
            "Aspergillus_niger_CBS_513.88|XP_001398067.1",
            "Aspergillus_awamori_IFM_58123|GCB24888.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_09150-RA",
            "Aspergillus_fumigatus_Af293|EAL84262.1",
            "Aspergillus_fischeri_NRRL181|XP_001267441.1",
            "Aspergillus_niger_CBS_513.88|XP_001400898.1",
            "Aspergillus_awamori_IFM_58123|GCB19337.1",
            "Aspergillus_niger_CBS_513.88|XP_001397349.2",
            "Aspergillus_awamori_IFM_58123|GCB25420.1",
            "Aspergillus_fumigatus_Af293|EAL87116.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_03570-RA",
            "Aspergillus_fischeri_NRRL181|XP_001265570.1",
            "Aspergillus_niger_CBS_513.88|XP_001390417.1",
            "Aspergillus_niger_CBS_513.88|XP_001390415.2",
            "Aspergillus_awamori_IFM_58123|GCB24160.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_08472-RA",
            "Aspergillus_fumigatus_Af293|EAL94045.1",
            "Aspergillus_fischeri_NRRL181|XP_001261225.1",
            "Aspergillus_niger_CBS_513.88|XP_001402298.1",
            "Aspergillus_awamori_IFM_58123|GCB23392.1",
            "Aspergillus_fumigatus_Af293|EAL93843.2",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_04905-RA",
            "Aspergillus_fischeri_NRRL181|XP_001261009.1",
            "Aspergillus_niger_CBS_513.88|XP_001397083.2",
            "Aspergillus_awamori_IFM_58123|GCB27653.1",
            "Aspergillus_niger_CBS_513.88|XP_001391581.1",
            "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate1",
            "Aspergillus_awamori_IFM_58123|GCB17486.1-duplicate",
            "Aspergillus_awamori_IFM_58123|GCB17486.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
            "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
            "Aspergillus_fumigatus_Af293|EAL85095.2",
            "Aspergillus_fischeri_NRRL181|XP_001261692.1",
            "Aspergillus_niger_CBS_513.88|XP_001401336.1",
            "Aspergillus_awamori_IFM_58123|GCB19008.1",
        ]

        prune_subtree(all_tips, terms, newtree)

        ## check results
        for term0, term1 in zip(tree.get_terminals(), newtree.get_terminals()):
            assert term0.name == term1.name
            assert term0.branch_length == term1.branch_length


class TestReadInputFiles(object):
    def test_read_input_files(self):
        ## setup
        expected_tree = Phylo.read(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            "newick",
        )
        expected_tree.root_at_midpoint()
        expected_fasta = SeqIO.to_dict(
            SeqIO.parse(
                f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
                "fasta",
            )
        )

        ## execution
        fasta = f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit"
        tree = (
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile"
        )

        (tree, fasta) = read_input_files(tree, fasta, False)

        ## check results
        for term0, term1 in zip(tree.get_terminals(), expected_tree.get_terminals()):
            assert term0.name == term1.name
            assert term0.branch_length == term1.branch_length
        for key, value in fasta.items():
            assert key in expected_fasta.keys()
            assert value.seq == expected_fasta[key].seq


# class TestWriteOutputFastaAndAccountForAssignedTipsSingleCopyCase(object):
#     def test_write_output_fasta_and_account_for_assigned_tips_single_copy_case(self):
#         ## setup
#         fasta = f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit"
#         fasta_dict = SeqIO.to_dict(
#             SeqIO.parse(
#                 f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
#                 "fasta",
#             )
#         )
#         subgroup_counter = 4
#         terms = [
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fischeri_NRRL181|XP_001261692.1",
#             "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#             "Aspergillus_awamori_IFM_58123|GCB19008.1",
#         ]
#         assigned_tips = [
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_06972-RA",
#             "Aspergillus_fumigatus_Af293|EAL85274.1",
#             "Aspergillus_fischeri_NRRL181|XP_001262055.1",
#             "Aspergillus_niger_CBS_513.88|XP_001398067.1",
#             "Aspergillus_awamori_IFM_58123|GCB24888.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_09150-RA",
#             "Aspergillus_fumigatus_Af293|EAL84262.1",
#             "Aspergillus_fischeri_NRRL181|XP_001267441.1",
#             "Aspergillus_niger_CBS_513.88|XP_001400898.1",
#             "Aspergillus_awamori_IFM_58123|GCB19337.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_08472-RA",
#             "Aspergillus_fumigatus_Af293|EAL94045.1",
#             "Aspergillus_fischeri_NRRL181|XP_001261225.1",
#             "Aspergillus_niger_CBS_513.88|XP_001402298.1",
#             "Aspergillus_awamori_IFM_58123|GCB23392.1",
#             "Aspergillus_fumigatus_Af293|EAL93843.2",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_04905-RA",
#             "Aspergillus_fischeri_NRRL181|XP_001261009.1",
#             "Aspergillus_niger_CBS_513.88|XP_001397083.2",
#             "Aspergillus_awamori_IFM_58123|GCB27653.1",
#         ]
#         expected_assigned_tips = [
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_06972-RA",
#             "Aspergillus_fumigatus_Af293|EAL85274.1",
#             "Aspergillus_fischeri_NRRL181|XP_001262055.1",
#             "Aspergillus_niger_CBS_513.88|XP_001398067.1",
#             "Aspergillus_awamori_IFM_58123|GCB24888.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_09150-RA",
#             "Aspergillus_fumigatus_Af293|EAL84262.1",
#             "Aspergillus_fischeri_NRRL181|XP_001267441.1",
#             "Aspergillus_niger_CBS_513.88|XP_001400898.1",
#             "Aspergillus_awamori_IFM_58123|GCB19337.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_08472-RA",
#             "Aspergillus_fumigatus_Af293|EAL94045.1",
#             "Aspergillus_fischeri_NRRL181|XP_001261225.1",
#             "Aspergillus_niger_CBS_513.88|XP_001402298.1",
#             "Aspergillus_awamori_IFM_58123|GCB23392.1",
#             "Aspergillus_fumigatus_Af293|EAL93843.2",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_04905-RA",
#             "Aspergillus_fischeri_NRRL181|XP_001261009.1",
#             "Aspergillus_niger_CBS_513.88|XP_001397083.2",
#             "Aspergillus_awamori_IFM_58123|GCB27653.1",
#             "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_05774-RA-duplicate",
#             "Aspergillus_fumigatus_Af293|EAL85095.2-duplicate",
#             "Aspergillus_fischeri_NRRL181|XP_001261692.1",
#             "Aspergillus_niger_CBS_513.88|XP_001401336.1",
#             "Aspergillus_awamori_IFM_58123|GCB19008.1",
#         ]

#         output_path = f"{here.parent.parent}/samples/"
#         expected_subgroup_counter = 5

#         place_holder = ""

#         inparalog_handling = dict()
#         inparalog_handling_summary = dict()

#         ## execution
#         (
#             subgroup_counter,
#             assigned_tips,
#         ) = write_output_fasta_and_account_for_assigned_tips_single_copy_case(
#             fasta,
#             subgroup_counter,
#             terms,
#             fasta_dict,
#             assigned_tips,
#             place_holder,
#             place_holder,
#             output_path,
#             inparalog_handling,
#             inparalog_handling_summary,
#         )

#         ## check results
#         assert subgroup_counter == expected_subgroup_counter
#         assert assigned_tips == expected_assigned_tips


class TestDetermineIfInputTreeIsSingleCopy(object):
    def test_determine_if_input_tree_is_single_copy_true(self):
        ## setup
        taxa = [
            "Aspergillus_fumigatus_Af293",
            "Aspergillus_oerlinghausenensis_CBS139183",
            "Aspergillus_fischeri_NRRL181",
        ]
        all_tips = [
            "Aspergillus_fumigatus_Af293|EAL84942.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_06972-RA",
            "Aspergillus_fischeri_NRRL181|XP_001262055.1",
        ]
        expected_res = True

        ## execution
        res = check_if_single_copy(taxa, all_tips)

        ## check results
        assert res == expected_res

    def test_determine_if_input_tree_is_single_copy_false(self):
        ## setup
        taxa = [
            "Aspergillus_fumigatus_Af293",
            "Aspergillus_oerlinghausenensis_CBS139183",
            "Aspergillus_fischeri_NRRL181",
        ]
        all_tips = [
            "Aspergillus_fumigatus_Af293|EAL84942.1",
            "Aspergillus_oerlinghausenensis_CBS139183|A_oerling_CBS139183_06972-RA",
            "Aspergillus_fischeri_NRRL181|XP_001262055.1",
            "Aspergillus_fischeri_NRRL181|XP_001262055.1",
        ]
        expected_res = False

        ## execution
        res = check_if_single_copy(taxa, all_tips)

        ## check results
        assert res == expected_res


class TestInparalogMedianSelection(object):
    def test_median_seq_len_more_than_two_duplicates(self):
        newtree = Phylo.read(StringIO("(sp|g1:0.1,sp|g2:0.2,sp|g3:0.3);"), "newick")
        dups = ["sp|g1", "sp|g2", "sp|g3"]
        terms = dups.copy()
        fasta_dict = {
            "sp|g1": SeqRecord(Seq("A"), id="sp|g1"),
            "sp|g2": SeqRecord(Seq("AA"), id="sp|g2"),
            "sp|g3": SeqRecord(Seq("AAA"), id="sp|g3"),
        }
        inparalog_handling = dict()

        _, updated_terms, inparalog_handling, pruned_tips = inparalog_to_keep_determination(
            newtree,
            fasta_dict,
            dups,
            terms,
            InparalogToKeep.median_seq_len,
            inparalog_handling,
        )

        assert updated_terms == ["sp|g2"]
        assert inparalog_handling["sp|g2"] == ["sp|g1", "sp|g3"]
        assert pruned_tips == ["sp|g1", "sp|g3"]

    def test_median_branch_len_more_than_two_duplicates(self):
        newtree = Phylo.read(StringIO("(sp|g1:0.1,sp|g2:0.2,sp|g3:0.3);"), "newick")
        dups = ["sp|g1", "sp|g2", "sp|g3"]
        terms = dups.copy()
        fasta_dict = {
            "sp|g1": SeqRecord(Seq("A"), id="sp|g1"),
            "sp|g2": SeqRecord(Seq("AA"), id="sp|g2"),
            "sp|g3": SeqRecord(Seq("AAA"), id="sp|g3"),
        }
        inparalog_handling = dict()

        _, updated_terms, inparalog_handling, pruned_tips = inparalog_to_keep_determination(
            newtree,
            fasta_dict,
            dups,
            terms,
            InparalogToKeep.median_branch_len,
            inparalog_handling,
        )

        assert updated_terms == ["sp|g2"]
        assert inparalog_handling["sp|g2"] == ["sp|g1", "sp|g3"]
        assert pruned_tips == ["sp|g1", "sp|g3"]
