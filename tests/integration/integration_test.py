import pytest
import sys

from pathlib import Path

from orthosnap.helper import InparalogToKeep
from orthosnap.orthosnap import execute


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
            snap_trees=False,
            inparalog_to_keep=InparalogToKeep.longest_seq_len,
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
