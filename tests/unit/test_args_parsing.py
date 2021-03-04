import pytest
from argparse import Namespace

from orthosnap.args_processing import process_args


@pytest.fixture
def args():
    kwargs = dict(
        tree="tests/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
        fasta="tests/samples/OG0000010.renamed.fa.mafft.clipkit",
        support=80,
        occupancy=1,
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

    def test_process_args_occupancy_above_range(self, args):
        args.occupancy = 1.1
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
