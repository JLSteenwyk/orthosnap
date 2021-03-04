import os
import pytest
import subprocess


class TestEntrypoint(object):
    def test_help(self):
        cmd = "orthosnap --help"
        exit_status = os.system(cmd)
        assert exit_status == 0

    def test_run(self):
        cmd = "orthosnap -t tests/samples/OG0000010.renamed.fa.mafft.clipkit.treefile -f tests/samples/OG0000010.renamed.fa.mafft.clipkit"
        exit_status = os.system(cmd)
        assert exit_status == 0

    def test_input_error0(self):
        cmd = "orthosnap -t tests/samples/OG0000010.renamed.fa.mafft.clipkit.treefile -f tests/samples/does_not_exist"
        response = subprocess.check_output([cmd], stderr=subprocess.STDOUT, shell=True)
        assert response == b"Input fasta does not exist\n"

    def test_input_error1(self):
        cmd = "orthosnap -t tests/samples/does_not_exist -f tests/samples/OG0000010.renamed.fa.mafft.clipkit"
        response = subprocess.check_output([cmd], stderr=subprocess.STDOUT, shell=True)
        assert response == b"Input tree does not exist\n"

    def test_run_no_args(self):
        cmd = "orthosnap"
        exit_status = os.system(cmd)
        assert exit_status == 0
