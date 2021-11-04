Usage
=====

**This section covers basic orthosnap usage. We recommend going through the 
tutorial as well, which gives a worked example of what orthosnap does.**

Orthosnap fits into a larger pipeline to construct phylogenomic data matrices.
After calling orthologous groups of genes, it is common for the number of 
single-copy orthologous groups of genes to be relatively few compared to all
orthologous groups. Typically, only single-copy orthologous groups of genes
are used for phylogenomic analyses. Orthosnap identifies subgroups of
single-copy orthologous genes in multi-copy orthologous groups of genes.

The input for orthosnap is the phylogeny of the orthologous group of genes
and the unaligned FASTA file of sequences used to infer the phylogeny. The
output of orthosnap are multi-FASTA files of subgroups of single-copy orthologous
genes appropriate for downstream analyses.

|

Basic usage
-----------

The following command is the simpliest iteration of orthosnap and will suffice
for most use cases:

.. code-block:: shell

   $ orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre

Here, two required arguments are shown, the -f/\-\-fasta argument, which specifies
the unaligned orthologous group of sequences, and the -t/\-\-tree argument, which
species the phylogeny inferred from the orthologous group of sequences. The fasta
file should be in FASTA format and the tree file should be in newick format.

|

Accounting for tree uncertainty
-------------------------------

As part of the orthosnap pipeline, in-paralogs and paralogs are pruned following the
approach described in `PhyloTreePruner <https://journals.sagepub.com/doi/10.4137/EBO.S12813>`_.
To do so, poorly collapsed bipartitions are collapsed to account for tree uncertainty.

The default threshold for collapsing bipartitions is 80 and can be modified using the 
-s/\-\-support argument. For example, bipartitions can be collapsed using a threshold
of 70 using the following command:

.. code-block:: shell

   $ orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre -s 70

|

All options
-----------

+-----------------------------+---------------------------------------------------------------------------------+
| Option                      | Usage and meaning                                                               |
+=============================+=================================================================================+
| -h/\-\-help                 | Print help message                                                              |
+-----------------------------+---------------------------------------------------------------------------------+
| -v/\-\-version              | Print software version                                                          |
+-----------------------------+---------------------------------------------------------------------------------+
| -t/\-\-tree                 | Input tree file (format: newick)                                                |
+-----------------------------+---------------------------------------------------------------------------------+
| -s/\-\-support              | Bipartition support threshold for collapsing uncertain branches (default: 80)   |
+-----------------------------+---------------------------------------------------------------------------------+
| -o/\-\-occupancy            | Occupancy threshold for identifying a subgroup of interest (default: 1.0)       |
+-----------------------------+---------------------------------------------------------------------------------+
*Note, although we provide a parameter for specify subgroup occupancy, we recommend not changing this parameter.
The -o/\-\-occupancy parameter was introduced into orthosnap for exploratory analyses by the developers.*