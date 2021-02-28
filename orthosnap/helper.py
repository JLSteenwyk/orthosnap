import sys

from Bio import Phylo
from Bio import SeqIO

def read_input_files(
    tree: str,
    fasta: str
):
    """
    create weights
    --------------
    to make each gene weigh the same (or contribute the same)
    the signal of each quartet from a given gene will be weighted
    """

    tree = Phylo.read(tree, "newick")

    # fasta = AlignIO.read(open(fasta), "fasta")
    fasta = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    return tree, fasta


