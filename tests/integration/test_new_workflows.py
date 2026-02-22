from pathlib import Path

import pytest

from orthosnap.orthosnap import main


HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
SAMPLE_TREE = ROOT / "samples" / "OG0000010.renamed.fa.mafft.clipkit.treefile"
SAMPLE_FASTA = ROOT / "samples" / "OG0000010.renamed.fa.mafft.clipkit"


@pytest.mark.integration
class TestNewWorkflows:
    def test_validate_only(self, tmp_path):
        main(
            [
                "-t",
                str(SAMPLE_TREE),
                "-f",
                str(SAMPLE_FASTA),
                "--validate-only",
                "-op",
                str(tmp_path),
            ]
        )

        assert not list(tmp_path.glob("*.orthosnap.*.fa"))

    def test_structured_output_and_resume(self, tmp_path):
        args = [
            "-t",
            str(SAMPLE_TREE),
            "-f",
            str(SAMPLE_FASTA),
            "--structured-output",
            "-op",
            str(tmp_path),
        ]

        main(args)

        run_json = tmp_path / f"{SAMPLE_FASTA.name}.orthosnap.run.json"
        subgroup_tsv = tmp_path / f"{SAMPLE_FASTA.name}.orthosnap.subgroups.tsv"
        first_subgroups = list(tmp_path.glob(f"{SAMPLE_FASTA.name}.orthosnap.*.fa"))

        assert run_json.exists()
        assert subgroup_tsv.exists()
        assert len(first_subgroups) > 0

        main(args + ["--resume"])

        second_subgroups = list(tmp_path.glob(f"{SAMPLE_FASTA.name}.orthosnap.*.fa"))
        assert len(second_subgroups) == len(first_subgroups)

    def test_manifest_mode(self, tmp_path):
        manifest = tmp_path / "manifest.tsv"
        manifest.write_text(
            "tree\tfasta\tid\n"
            f"{SAMPLE_TREE}\t{SAMPLE_FASTA}\tjob1\n"
        )

        out_dir = tmp_path / "batch"
        main([
            "--manifest",
            str(manifest),
            "--structured-output",
            "-op",
            str(out_dir),
        ])

        summary_tsv = list(out_dir.glob("manifest_summary_*.tsv"))
        summary_json = list(out_dir.glob("manifest_summary_*.json"))

        assert len(summary_tsv) == 1
        assert len(summary_json) == 1
        assert (out_dir / "job1").exists()

    def test_bootstrap_consensus_mode(self, tmp_path):
        bootstrap = tmp_path / "bootstrap_trees.txt"
        bootstrap.write_text(f"{SAMPLE_TREE}\n{SAMPLE_TREE}\n")

        main(
            [
                "-t",
                str(SAMPLE_TREE),
                "-f",
                str(SAMPLE_FASTA),
                "--bootstrap-trees",
                str(bootstrap),
                "--consensus-min-frequency",
                "0.5",
                "--structured-output",
                "-op",
                str(tmp_path),
            ]
        )

        consensus_tsv = tmp_path / f"{SAMPLE_FASTA.name}.orthosnap.consensus.tsv"
        consensus_fa = list(tmp_path.glob(f"{SAMPLE_FASTA.name}.orthosnap.consensus_*.fa"))

        assert consensus_tsv.exists()
        assert len(consensus_fa) > 0

    def test_occupancy_fraction_and_count_modes(self, tmp_path):
        out_fraction = tmp_path / "fraction"
        out_count = tmp_path / "count"

        main(
            [
                "-t",
                str(SAMPLE_TREE),
                "-f",
                str(SAMPLE_FASTA),
                "--occupancy-fraction",
                "0.5",
                "--structured-output",
                "-op",
                str(out_fraction),
            ]
        )

        main(
            [
                "-t",
                str(SAMPLE_TREE),
                "-f",
                str(SAMPLE_FASTA),
                "--occupancy-count",
                "3",
                "--structured-output",
                "-op",
                str(out_count),
            ]
        )

        assert (out_fraction / f"{SAMPLE_FASTA.name}.orthosnap.run.json").exists()
        assert (out_count / f"{SAMPLE_FASTA.name}.orthosnap.run.json").exists()
