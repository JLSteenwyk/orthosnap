import pytest
import sys

from pathlib import Path

from orthosnap.orthosnap import execute

here = Path(__file__)


@pytest.mark.integration
class TestIntegration(object):
    def test_simple_no_change(self):
        """
        """
        kwargs = dict(
            tree=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.treefile",
            fasta=f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit",
            support=80,
            occupancy=1,
        )
        execute(**kwargs)

        with open(
            f"{here.parent.parent}/expected/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa", "r"
        ) as expected:
            expected_content = expected.read()

        with open(f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.0.fa", "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa", "r"
        ) as expected:
            expected_content = expected.read()

        with open(f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.1.fa", "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa", "r"
        ) as expected:
            expected_content = expected.read()

        with open(f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.2.fa", "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa", "r"
        ) as expected:
            expected_content = expected.read()

        with open(f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.3.fa", "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content

        with open(
            f"{here.parent.parent}/expected/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa", "r"
        ) as expected:
            expected_content = expected.read()

        with open(f"{here.parent.parent}/samples/OG0000010.renamed.fa.mafft.clipkit.orthosnap.4.fa", "r") as out_file:
            output_content = out_file.read()

        assert expected_content == output_content
