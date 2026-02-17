import pytest
from pathlib import Path

from orthosnap.helper import InparalogToKeep
from orthosnap.orthosnap import execute, main


here = Path(__file__)


@pytest.mark.integration
class TestIntegration(object):
    def test_default_param_OG0000010(self):
        """"""
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=80,
            occupancy=5,
            rooted=False,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{here.parent.parent}/samples/",
            report_inparalog_handling=False,
            delimiter="|",
        )
        execute(**kwargs)

        with open(
            f"{here.parent.parent}/expected/test_default_param_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_default_param_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_default_param_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_default_param_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_default_param_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_occupancy_two_tenths_OG0000010(self):
        """"""
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=80,
            occupancy=0.2,
            rooted=False,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{here.parent.parent}/samples/",
            report_inparalog_handling=False,
            delimiter="|",
        )
        execute(**kwargs)

        with open(
            f"{here.parent.parent}/expected/test_occupancy_two_tenths_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_occupancy_two_tenths_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_occupancy_two_tenths_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_occupancy_two_tenths_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_occupancy_two_tenths_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_occupancy_two_tenths_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.5.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.5.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_occupancy_two_tenths_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.6.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.6.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_occupancy_two_tenths_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.7.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.7.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_occupancy_two_tenths_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.8.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.8.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_occupancy_two_tenths_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.9.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.9.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_support_value_60_OG0000010(self):
        """"""
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=60,
            occupancy=5,
            rooted=False,
            snap_trees=True,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{here.parent.parent}/samples/",
            report_inparalog_handling=False,
            delimiter="|",
        )
        execute(**kwargs)

        with open(
            f"{here.parent.parent}/expected/test_support_value_60_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_support_value_60_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_support_value_60_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_support_value_60_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_support_value_60_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_snap_trees_argument(self):
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=80,
            occupancy=5,
            rooted=False,
            snap_trees=True,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{here.parent.parent}/samples/",
            report_inparalog_handling=False,
            delimiter="|",
        )
        execute(**kwargs)

        with open(
            f"{here.parent.parent}/expected/test_snap_trees_argument/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_snap_trees_argument/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_snap_trees_argument/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_snap_trees_argument/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_snap_trees_argument/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_snap_trees_argument/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.tre",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.tre",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_snap_trees_argument/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.tre",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.tre",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_snap_trees_argument/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.tre",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.tre",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_snap_trees_argument/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.tre",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.tre",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_snap_trees_argument/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.tre",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.tre",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_rooted_argument(self):
        """"""
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=60,
            occupancy=5,
            rooted=True,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{here.parent.parent}/samples/",
            report_inparalog_handling=False,
            delimiter="|",
        )
        execute(**kwargs)

        with open(
            f"{here.parent.parent}/expected/test_rooted_argument/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_rooted_argument/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_rooted_argument/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_rooted_argument/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_rooted_argument/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_already_single_copy(self):
        """"""
        testargs = [
            "orthosnap",
            "-t",
            f"{here.parent.parent}/samples/already_single_copy.tre",
            "-f",
            f"{here.parent.parent}/samples/already_single_copy.fa",
            "-o",
            "5",
        ]

        with pytest.raises(SystemExit) as pytest_wrapped_e:
            main()

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

    def test_shortest_branch_len(self):
        """"""
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=80,
            occupancy=3,
            rooted=False,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.shortest_branch_len,
            output_path=f"{here.parent.parent}/samples/",
            report_inparalog_handling=False,
            delimiter="|",
        )
        execute(**kwargs)

        with open(
            f"{here.parent.parent}/expected/test_shortest_branch_len/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_shortest_branch_len/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_shortest_branch_len/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_shortest_branch_len/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_shortest_branch_len/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_median_branch_len(self):
        """"""
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=80,
            occupancy=3,
            rooted=False,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.shortest_branch_len,
            output_path=f"{here.parent.parent}/samples/",
            report_inparalog_handling=False,
            delimiter="|",
        )
        execute(**kwargs)

        with open(
            f"{here.parent.parent}/expected/test_median_branch_len/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_median_branch_len/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_median_branch_len/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_median_branch_len/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_median_branch_len/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_specifying_output_directory(self):
        """"""
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=80,
            occupancy=5,
            rooted=False,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{here.parent.parent}/samples/specified_dir/",
            report_inparalog_handling=False,
            delimiter="|",
        )
        execute(**kwargs)

        with open(
            f"{here.parent.parent}/expected/test_default_param_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/specified_dir/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_default_param_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/specified_dir/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_default_param_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/specified_dir/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_default_param_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/specified_dir/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_default_param_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/specified_dir/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_specifying_output_directory_without_ending_slash(self):
        """"""
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=80,
            occupancy=5,
            rooted=False,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{here.parent.parent}/samples/specified_dir",
            report_inparalog_handling=False,
            delimiter="|",
        )
        execute(**kwargs)

        with open(
            f"{here.parent.parent}/expected/test_default_param_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/specified_dir/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_default_param_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/specified_dir/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_default_param_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/specified_dir/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_default_param_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/specified_dir/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_default_param_OG0000010/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/specified_dir/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_inparalog_summary_file(self):
        """"""
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=80,
            occupancy=3,
            rooted=False,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{here.parent.parent}/samples/",
            report_inparalog_handling=True,
            delimiter="|",
        )

        execute(**kwargs)

        with open(
            f"{here.parent.parent}/expected/OG0000010.renamed.fa.mafft.clipkit.inparalog_report.txt",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.inparalog_report.txt",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_default_param_fake_ogs_inparalog_report_occupancy_1(self):
        """"""
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/dataset/fake_orthologous_group_of_genes_tree.tre",
            fasta=f"{here.parent.parent}/samples/dataset/fake_orthologous_group_of_genes.faa",
            support=80,
            occupancy=1.0,
            rooted=False,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{here.parent.parent}/samples/dataset/",
            report_inparalog_handling=False,
            delimiter="|",
        )
        execute(**kwargs)

        with open(
            f"{here.parent.parent}/expected/fake_ogs/fake_orthologous_group_of_genes.faa.orthosnap.0.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/dataset/fake_orthologous_group_of_genes.faa.orthosnap.0.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/fake_ogs/fake_orthologous_group_of_genes.faa.orthosnap.1.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/dataset/fake_orthologous_group_of_genes.faa.orthosnap.1.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/fake_ogs/fake_orthologous_group_of_genes.faa.orthosnap.2.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/dataset/fake_orthologous_group_of_genes.faa.orthosnap.2.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/fake_ogs/fake_orthologous_group_of_genes.faa.orthosnap.3.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/dataset/fake_orthologous_group_of_genes.faa.orthosnap.3.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/fake_ogs/fake_orthologous_group_of_genes.faa.orthosnap.4.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/dataset/fake_orthologous_group_of_genes.faa.orthosnap.4.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_alt_delimiter_OG0000005(self):
        """"""
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000005.ampersand_delimiter.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000005.ampersand_delimiter.fa",
            support=80,
            occupancy=5,
            rooted=False,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{here.parent.parent}/samples/",
            report_inparalog_handling=False,
            delimiter="&",
        )
        execute(**kwargs)

        with open(
            f"{here.parent.parent}/expected/test_alt_delimiter_OG0000005/OG0000005.ampersand_delimiter.fa.orthosnap.0.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000005.ampersand_delimiter.fa.orthosnap.0.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/test_alt_delimiter_OG0000005/OG0000005.ampersand_delimiter.fa.orthosnap.1.fa",
            "r",
        ) as expected:
            expected_content = expected.read()

        with open(
            f"{here.parent.parent}/samples/OG0000005.ampersand_delimiter.fa.orthosnap.1.fa",
            "r",
        ) as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

    def test_no_inparalog_report_written_when_disabled(self, tmp_path):
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=80,
            occupancy=5,
            rooted=False,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{tmp_path}/",
            report_inparalog_handling=False,
            delimiter="|",
        )
        execute(**kwargs)

        report_file = tmp_path / "OG0000010.renamed.fa.mafft.clipkit.inparalog_report.txt"
        output_fasta = tmp_path / "OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa"

        assert output_fasta.exists()
        assert not report_file.exists()

    def test_execute_output_path_without_trailing_slash(self, tmp_path):
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=80,
            occupancy=5,
            rooted=False,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{tmp_path}",
            report_inparalog_handling=True,
            delimiter="|",
        )
        execute(**kwargs)

        report_file = tmp_path / "OG0000010.renamed.fa.mafft.clipkit.inparalog_report.txt"
        output_fasta = tmp_path / "OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa"

        assert output_fasta.exists()
        assert report_file.exists()

    def test_occupancy_boundary_behavior(self, tmp_path):
        eq_path = tmp_path / "occ_eq"
        gt_path = tmp_path / "occ_gt"
        eq_path.mkdir()
        gt_path.mkdir()

        eq_kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=80,
            occupancy=5,
            rooted=False,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{eq_path}",
            report_inparalog_handling=False,
            delimiter="|",
        )
        execute(**eq_kwargs)
        eq_count = len(list(eq_path.glob("*.orthosnap.*.fa")))
        assert eq_count == 5

        gt_kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=80,
            occupancy=5.01,
            rooted=False,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{gt_path}",
            report_inparalog_handling=False,
            delimiter="|",
        )
        execute(**gt_kwargs)
        gt_count = len(list(gt_path.glob("*.orthosnap.*.fa")))
        assert gt_count == 0

    def test_mixed_delimiter_input_raises_clean_error(self, tmp_path, capsys):
        tree_path = tmp_path / "mixed_delim.tree"
        fasta_path = tmp_path / "mixed_delim.fa"
        tree_path.write_text("(sp1|g1:0.1,sp2g2:0.1);\n")
        fasta_path.write_text(">sp1|g1\nAAAA\n>sp2g2\nAAAA\n")

        kwargs = dict(
            tree=f"{tree_path}",
            fasta=f"{fasta_path}",
            support=80,
            occupancy=1,
            rooted=True,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{tmp_path}",
            report_inparalog_handling=False,
            delimiter="|",
        )

        with pytest.raises(SystemExit):
            execute(**kwargs)

        std_out = capsys.readouterr().out
        assert "ERROR: Delimiter does not exist in FASTA headers." in std_out

    def test_plot_snap_ogs_output_file_created(self, tmp_path):
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=80,
            occupancy=5,
            rooted=False,
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
            output_path=f"{tmp_path}",
            report_inparalog_handling=False,
            delimiter="|",
            plot_snap_ogs_output=True,
            plot_format="png",
        )
        execute(**kwargs)

        plot_file = tmp_path / "OG0000010.renamed.fa.mafft.clipkit.orthosnap.subgroups.png"
        assert plot_file.exists()
