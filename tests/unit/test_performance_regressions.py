from collections import Counter
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

import pytest
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from orthosnap.helper import (
    InparalogToKeep,
    build_subtree_taxa_cache,
    clone_induced_tree,
    handle_multi_copy_subtree,
    handle_single_copy_subtree,
)
from orthosnap.orthosnap import _extract_subgroups, main


SAMPLES = Path(__file__).resolve().parents[1] / "samples"
SAMPLE_TREE = SAMPLES / "OG0000010.renamed.fa.mafft.clipkit.treefile"
SAMPLE_FASTA = SAMPLES / "OG0000010.renamed.fa.mafft.clipkit"


def _extract_with_legacy_multi_copy(tree, fasta_dict, mode, support):
    assigned_tips = set()
    subgroup_counter = 0
    subgroup_records = []
    inparalog_handling = {}
    inparalog_handling_summary = {}
    cache = build_subtree_taxa_cache(tree, "|")

    for clade in cache.internal_clades[1:]:
        terms, terms_set, taxa_counts, counts = cache[clade]
        if len(taxa_counts) < 5 or not assigned_tips.isdisjoint(terms_set):
            continue

        if set(counts) == {1}:
            (
                subgroup_counter,
                assigned_tips,
                inparalog_handling,
                inparalog_handling_summary,
            ) = handle_single_copy_subtree(
                clade,
                terms,
                subgroup_counter,
                "synthetic.fasta",
                support,
                fasta_dict,
                assigned_tips,
                False,
                "",
                inparalog_handling,
                inparalog_handling_summary,
                False,
                subgroup_records,
                False,
            )
        else:
            (
                subgroup_counter,
                assigned_tips,
                inparalog_handling,
                inparalog_handling_summary,
            ) = handle_multi_copy_subtree(
                clade,
                terms,
                subgroup_counter,
                "synthetic.fasta",
                support,
                fasta_dict,
                assigned_tips,
                taxa_counts,
                False,
                mode,
                "",
                inparalog_handling,
                inparalog_handling_summary,
                False,
                "|",
                subgroup_records,
                False,
            )

    return subgroup_counter, subgroup_records


@pytest.mark.parametrize("support", [60, 80])
@pytest.mark.parametrize("mode", list(InparalogToKeep))
def test_indexed_extraction_matches_legacy_path(mode, support):
    tree = Phylo.read(SAMPLE_TREE, "newick")
    tree.root_at_midpoint()
    fasta_dict = SeqIO.to_dict(SeqIO.parse(SAMPLE_FASTA, "fasta"))

    expected = _extract_with_legacy_multi_copy(
        deepcopy(tree),
        fasta_dict,
        mode,
        support,
    )
    result = _extract_subgroups(
        tree=tree,
        fasta="synthetic.fasta",
        fasta_dict=fasta_dict,
        support=support,
        occupancy=5,
        snap_trees=False,
        inparalog_to_keep=mode,
        output_path="",
        report_inparalog_handling=False,
        delimiter="|",
        write_outputs=False,
    )

    assert (result["subgroup_counter"], result["subgroup_records"]) == expected


def _balanced_tree(tip_count):
    clades = [
        Clade(name=f"taxon_{index}|gene_{index}", branch_length=0.1)
        for index in range(tip_count)
    ]
    while len(clades) > 1:
        clades = [
            Clade(clades=clades[index : index + 2], confidence=100)
            for index in range(0, len(clades), 2)
        ]
    return Tree(root=clades[0])


def _unbalanced_tree(tip_count):
    root = Clade(name="taxon_0|gene_0", branch_length=0.1)
    for index in range(1, tip_count):
        root = Clade(
            confidence=100,
            clades=[
                root,
                Clade(name=f"taxon_{index}|gene_{index}", branch_length=0.1),
            ],
        )
    return Tree(root=root)


def _multi_group_tree(group_count, taxa_per_group, shape):
    subgroup_clades = []
    fasta_dict = {}
    expected_groups = []
    for group_index in range(group_count):
        names = [
            f"taxon_{taxon_index}|gene_{group_index}"
            for taxon_index in range(taxa_per_group)
        ]
        expected_groups.append(names)
        tips = []
        for name in names:
            tips.append(Clade(name=name, branch_length=0.1))
            fasta_dict[name] = SeqRecord(Seq("MAAAAA"), id=name, description="")
        subgroup_clades.append(Clade(confidence=100, clades=tips))

    if shape == "balanced":
        current = subgroup_clades
        while len(current) > 1:
            current = [
                Clade(confidence=100, clades=current[index : index + 2])
                for index in range(0, len(current), 2)
            ]
        root = current[0]
    else:
        root = subgroup_clades[0]
        for subgroup in subgroup_clades[1:]:
            root = Clade(confidence=100, clades=[root, subgroup])
    return Tree(root=root), fasta_dict, expected_groups


@pytest.mark.parametrize("tree_factory", [_balanced_tree, _unbalanced_tree])
def test_compact_cache_preserves_order_and_linear_metadata(tree_factory):
    tree = tree_factory(256)
    cache = build_subtree_taxa_cache(tree, "|")

    assert cache.terminal_names == [tip.name for tip in tree.get_terminals()]
    assert len(cache.internal_clades) == len(tree.get_nonterminals())
    assert len(cache._summaries) == len(cache.internal_clades)
    assert not any(isinstance(value, Counter) for value in cache._summaries.values())

    for clade in cache.internal_clades[:: max(1, len(cache.internal_clades) // 8)]:
        terms, terms_set, taxa_counts, counts = cache[clade]
        expected_terms = [tip.name for tip in clade.get_terminals()]
        assert terms == expected_terms
        assert terms_set == set(expected_terms)
        assert sum(taxa_counts.values()) == len(expected_terms)
        assert counts == list(taxa_counts.values())


def test_compact_cache_handles_tree_deeper_than_python_recursion_limit():
    tree = _unbalanced_tree(1500)
    cache = build_subtree_taxa_cache(tree, "|")

    assert len(cache.terminal_names) == 1500
    assert len(cache.internal_clades) == 1499


@pytest.mark.parametrize("shape", ["balanced", "unbalanced"])
def test_indexed_extraction_on_synthetic_tree_shapes(shape):
    tree, fasta_dict, expected_groups = _multi_group_tree(24, 6, shape)
    result = _extract_subgroups(
        tree=tree,
        fasta="synthetic.fasta",
        fasta_dict=fasta_dict,
        support=80,
        occupancy=6,
        snap_trees=False,
        inparalog_to_keep=InparalogToKeep.longest_seq_len,
        output_path="",
        report_inparalog_handling=False,
        delimiter="|",
        write_outputs=False,
    )

    assert result["subgroup_counter"] == 24
    assert [record["tips"] for record in result["subgroup_records"]] == expected_groups


def test_cli_reuses_tree_and_fasta_parsed_during_validation(tmp_path):
    original_tree_read = Phylo.read
    original_fasta_parse = SeqIO.parse
    with patch(
        "orthosnap.orthosnap.Phylo.read",
        side_effect=original_tree_read,
    ) as tree_read, patch(
        "orthosnap.orthosnap.SeqIO.parse",
        side_effect=original_fasta_parse,
    ) as fasta_parse:
        main(
            [
                "-t",
                str(SAMPLE_TREE),
                "-f",
                str(SAMPLE_FASTA),
                "-op",
                str(tmp_path),
            ]
        )

    assert tree_read.call_count == 1
    assert fasta_parse.call_count == 1


def test_bootstrap_reuses_fasta_and_parses_each_tree_once(tmp_path):
    bootstrap_path = tmp_path / "bootstrap_trees.txt"
    bootstrap_path.write_text(f"{SAMPLE_TREE}\n{SAMPLE_TREE}\n")
    original_tree_read = Phylo.read
    original_fasta_parse = SeqIO.parse
    with patch(
        "orthosnap.orthosnap.Phylo.read",
        side_effect=original_tree_read,
    ) as tree_read, patch(
        "orthosnap.orthosnap.SeqIO.parse",
        side_effect=original_fasta_parse,
    ) as fasta_parse:
        main(
            [
                "-t",
                str(SAMPLE_TREE),
                "-f",
                str(SAMPLE_FASTA),
                "-op",
                str(tmp_path / "output"),
                "--bootstrap-trees",
                str(bootstrap_path),
            ]
        )

    assert tree_read.call_count == 3
    assert fasta_parse.call_count == 1


@pytest.mark.parametrize("selection", ["prefix", "alternating", "singleton"])
def test_induced_clone_matches_sequential_pruning(selection):
    tree = Phylo.read(SAMPLE_TREE, "newick")
    tree.root_at_midpoint()
    terminal_names = [tip.name for tip in tree.get_terminals()]
    if selection == "prefix":
        keep_tips = set(terminal_names[:8])
    elif selection == "alternating":
        keep_tips = set(terminal_names[::3])
    else:
        keep_tips = {terminal_names[len(terminal_names) // 2]}

    sequential = deepcopy(tree)
    for terminal in list(sequential.get_terminals()):
        if terminal.name not in keep_tips:
            sequential.prune(terminal)

    induced = clone_induced_tree(tree, keep_tips)
    assert induced.format("newick") == sequential.format("newick")
