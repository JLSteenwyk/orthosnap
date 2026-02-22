import pytest

from orthosnap.parser import create_parser


@pytest.fixture
def parser():
    return create_parser()


class TestParser(object):
    def test_required_only(self, parser):
        fasta = "my/input/file.fa"
        tree = "my/input/tree.tree"

        parsed = parser.parse_args(["-f", fasta, "-t", tree])

        assert parsed.fasta == fasta
        assert parsed.tree == tree

    def test_occupancy(self, parser):
        fasta = "my/input/file.fa"
        tree = "my/input/tree.tree"
        occupancy = ".2"

        parsed = parser.parse_args(["-f", fasta, "-t", tree, "-o", occupancy])

        assert parsed.fasta == fasta
        assert parsed.tree == tree
        assert parsed.occupancy == float(occupancy)

    def test_support(self, parser):
        fasta = "my/input/file.fa"
        tree = "my/input/tree.tree"
        support = "70"

        parsed = parser.parse_args(["-f", fasta, "-t", tree, "-s", support])

        assert parsed.fasta == fasta
        assert parsed.tree == tree
        assert parsed.support == float(support)

    def test_parser_when_no_args(self, parser):
        parsed = parser.parse_args([])
        assert parsed.tree is None
        assert parsed.fasta is None

    def test_plot_flags(self, parser):
        parsed = parser.parse_args(
            [
                "-f",
                "my/input/file.fa",
                "-t",
                "my/input/tree.tree",
                "-ps",
                "-pf",
                "pdf",
            ]
        )
        assert parsed.plot_snap_ogs is True
        assert parsed.plot_format == "pdf"

    def test_new_workflow_flags(self, parser):
        parsed = parser.parse_args(
            [
                "--manifest",
                "manifest.tsv",
                "--validate-only",
                "--resume",
                "--structured-output",
                "--bootstrap-trees",
                "bootstrap_paths.txt",
                "--consensus-min-frequency",
                "0.7",
                "--occupancy-count",
                "4",
            ]
        )
        assert parsed.manifest == "manifest.tsv"
        assert parsed.validate_only is True
        assert parsed.resume is True
        assert parsed.structured_output is True
        assert parsed.bootstrap_trees == "bootstrap_paths.txt"
        assert parsed.consensus_min_frequency == 0.7
        assert parsed.occupancy_count == 4
