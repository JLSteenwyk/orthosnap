import pytest

from argparse import Namespace

from orthosnap.args_processing import (
    count_unique_taxa_in_fasta,
    determine_occupancy_threshold,
    process_args,
    proper_round,
)
from orthosnap.helper import InparalogToKeep


@pytest.fixture
def args():
    kwargs = dict(
        tree="tests/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
        fasta="tests/samples/OG0000010.renamed.fa.mafft.clipkit",
        support=80,
        occupancy=1,
        rooted=False,
        snap_trees=False,
        inparalog_to_keep=InparalogToKeep.longest_seq_len,
        output_path="./tests/samples/",
        report_inparalog_handling=False,
        delimiter="|",
        plot_snap_ogs=False,
        plot_format="png",
        manifest=None,
        validate_only=False,
        resume=False,
        structured_output=False,
        bootstrap_trees=None,
        consensus_min_frequency=0.5,
        occupancy_count=None,
        occupancy_fraction=None,
    )
    return Namespace(**kwargs)


class TestArgsProcessing(object):
    def test_process_args_tree_file_dne(self, args):
        args.tree = "some/file/that/doesnt/exist"
        with pytest.raises(SystemExit):
            process_args(args)

    def test_process_args_fasta_file_dne(self, args):
        args.fasta = "some/file/that/doesnt/exist"
        with pytest.raises(SystemExit):
            process_args(args)

    def test_process_args_occupancy_below_range(self, args):
        args.occupancy = -0.1
        with pytest.raises(SystemExit):
            process_args(args)

    def test_process_args_support_above_range(self, args):
        args.support = 101
        with pytest.raises(SystemExit):
            process_args(args)

    def test_process_args_support_below_range(self, args):
        args.support = -0.1
        with pytest.raises(SystemExit):
            process_args(args)

    def test_determine_occupancy_threshold(self, args):
        res = determine_occupancy_threshold(args.fasta, args.delimiter)
        assert res == 3

    def test_count_unique_taxa_in_fasta(self, args):
        res = count_unique_taxa_in_fasta(args.fasta, args.delimiter)
        assert res == 5

    def test_proper_round0(self, args):
        res = proper_round(2.5)
        assert res == 3.0

    def test_proper_round1(self, args):
        res = proper_round(2.2)
        assert res == 2.0

    def test_rooted_arg(self, args):
        res = process_args(args)
        assert res["rooted"] == False

    def test_rooted_arg_true(self, args):
        args.rooted = True
        res = process_args(args)
        assert res["rooted"] == True

    def test_snap_trees(self, args):
        assert args.snap_trees == False

    def test_snap_trees_true(self, args):
        args.snap_trees = True
        assert args.snap_trees == True

    def test_inparalog_to_keep_default(self, args):
        assert args.inparalog_to_keep == InparalogToKeep.longest_seq_len

    def test_inparalog_to_keep_none(self, args):
        args.inparalog_to_keep = None
        res = process_args(args)
        assert res["inparalog_to_keep"] == InparalogToKeep.longest_seq_len

    def test_inparalog_to_keep_median_seq_len(self, args):
        res = InparalogToKeep.median_seq_len
        assert res == InparalogToKeep.median_seq_len

    def test_inparalog_to_keep_shortest_seq_len(self, args):
        res = InparalogToKeep.shortest_seq_len
        assert res == InparalogToKeep.shortest_seq_len

    def test_inparalog_to_keep_shortest_branch_len(self, args):
        res = InparalogToKeep.shortest_branch_len
        assert res == InparalogToKeep.shortest_branch_len

    def test_inparalog_to_keep_median_branch_len(self, args):
        res = InparalogToKeep.median_branch_len
        assert res == InparalogToKeep.median_branch_len

    def test_inparalog_to_keep_longest_branch_len(self, args):
        res = InparalogToKeep.longest_branch_len
        assert res == InparalogToKeep.longest_branch_len

    def test_output_path(self, args):
        args.output_path = "./tests/samples/"
        res = process_args(args)
        assert res["output_path"] == "./tests/samples/"

    def test_output_path_none(self, args):
        args.output_path = None
        args.fasta = "tests/expected/test_support_value_60_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa"
        res = process_args(args)
        assert res["output_path"] == "tests/expected/test_support_value_60_OG0000010/"

    def test_output_path_no_slash(self, args):
        args.output_path = "./tests/samples"
        res = process_args(args)
        assert res["output_path"] == "./tests/samples/"

    def test_output_path_none(self, args):
        args.fasta = "requirements.txt"  # fake stand in file
        args.output_path = None
        res = process_args(args)
        assert res["output_path"] == "./"

    def test_delimiter_default(self, args):
        res = process_args(args)
        assert res["delimiter"] == "|"

    def test_delimiter_custom(self, args):
        args.delimiter = "&"
        res = process_args(args)
        assert res["delimiter"] == "&"

    def test_plot_snap_ogs_default(self, args):
        res = process_args(args)
        assert res["plot_snap_ogs_output"] == False

    def test_plot_snap_ogs_true(self, args):
        args.plot_snap_ogs = True
        res = process_args(args)
        assert res["plot_snap_ogs_output"] == True

    def test_plot_format_default(self, args):
        args.plot_format = None
        res = process_args(args)
        assert res["plot_format"] == "png"

    def test_plot_format_custom(self, args):
        args.plot_format = "svg"
        res = process_args(args)
        assert res["plot_format"] == "svg"

    def test_occupancy_count_mode(self, args):
        args.occupancy = None
        args.occupancy_count = 4
        res = process_args(args)
        assert res["occupancy_mode"] == "count"
        assert res["occupancy"] == 4

    def test_occupancy_fraction_mode(self, args):
        args.occupancy = None
        args.occupancy_fraction = 0.5
        res = process_args(args)
        assert res["occupancy_mode"] == "fraction"
        assert res["occupancy"] == 3

    def test_occupancy_modes_conflict(self, args):
        args.occupancy = 3
        args.occupancy_count = 4
        with pytest.raises(SystemExit):
            process_args(args)

    def test_manifest_mode_without_tree_fasta(self, args, tmp_path):
        manifest = tmp_path / "manifest.tsv"
        manifest.write_text("tree\tfasta\nx\ty\n")
        args.tree = None
        args.fasta = None
        args.manifest = str(manifest)
        res = process_args(args)
        assert res["manifest"] == str(manifest)
