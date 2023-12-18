Usage
=====

**This section covers basic OrthoSNAP usage. We recommend going through the 
tutorial as well, which gives a worked example of what OrthoSNAP does.**

OrthoSNAP fits into a larger pipeline to construct phylogenomic data matrices.
After calling orthologous groups of genes, it is common for the number of 
single-copy orthologous groups of genes to be relatively few compared to all
orthologous groups. Typically, only single-copy orthologous groups of genes
are used for phylogenomic analyses. OrthoSNAP identifies subgroups of
single-copy orthologous genes in multi-copy orthologous groups of genes.

The input for OrthoSNAP is the phylogeny of the orthologous group of genes
and the unaligned FASTA file of sequences used to infer the phylogeny. The
output of OrthoSNAP are multi-FASTA files of subgroups of single-copy orthologous
genes appropriate for downstream analyses.

|

Basic usage
-----------

The following command is the simpliest iteration of OrthoSNAP and will suffice
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

As part of the OrthoSNAP pipeline, species-specific inparalogs and paralogs are pruned following the
approach described in `PhyloTreePruner <https://journals.sagepub.com/doi/10.4137/EBO.S12813>`_.
To do so, poorly collapsed bipartitions are collapsed to account for tree uncertainty.

The default threshold for collapsing bipartitions is 80 and can be modified using the 
-s/\-\-support argument. For example, bipartitions can be collapsed using a threshold
of 70 using the following command:

.. code-block:: shell

   $ orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre -s 70

|

Specifying which inparalog to keep
----------------------------------

During species-specific inparalog pruning, it is standard practice in transcriptomic analysis to keep the longest inparalog. However, this
may not always be the user's preference. OrthoSNAP allows flexibility in which inparalog should be kept using the -ip /\-\-inparalog_to_keep parameter.
Specifically, user's can choose to keep the longest sequence length (default; longest_seq_len), shortest sequence length (shortest_seq_len), or
median sequence length (median_seq_len) in the case of three or more inparalogs. User's can also choose which inparalog to 
keep based on tree-based metrics: the longest branch length (longest_branch_len), the shortest branch length (shortest_branch_len), or the
median branch length (median_branch_len) in the case of three or more inparalogs.

For example, the inparalog with the shortest branch length can be kept using the following command:

.. code-block:: shell

   $ orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre -ip shortest_branch_len

|

or the inparalog with the median sequence length can be kept using the following command:

.. code-block:: shell

   $ orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre -ip median_seq_len

|

Again, following transcriptomics, the default is to keep the longest sequence because (at least in theory)
it is the most complete gene annotation.

Report inparalog handling
-------------------------
To report inparalogs and specify which was kept per SNAP-OG, use the -rih, \-\-report_inparalog_handling
argument. The resulting file, which will have the suffix ".inparalog_report.txt," will have three columns: |br|
- col 1 is the orthogroup file |br|
- col 2 is the inparalog that was kept |br|
- col 3 is/are the inparalog/s that were trimmed separated by a semi-colon ";" |br|

To generate this file, use the following command:

.. code-block:: shell

   $ orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre -rih

|

All options
-----------

+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| Option                              | Usage and meaning                                                                                                                            |
+=====================================+==============================================================================================================================================+
| -h/\-\-help                         | Print help message                                                                                                                           |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| -v/\-\-version                      | Print software version                                                                                                                       |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| -t/\-\-tree                         | Input tree file (format: newick)                                                                                                             |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| -s/\-\-support                      | Bipartition support threshold for collapsing uncertain branches (default: 80)                                                                |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| -o/\-\-occupancy                    | Occupancy threshold for identifying a subgroup of interest (default: 50%)                                                                    |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| -r/\-\-roooted                      | boolean argument for whether the input phylogeny is already rooted (default: false)                                                          |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| -st/\-\-snap_trees                  | boolean argument for whether trees of SNAP-OGs should be outputted (default: false)                                                          |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| -ip/\-\-inparalog_to_keep           | determine which sequence to keep in the case of species-specific inparalogs using sequence- or tree-based options (default: longest_seq_len) |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| -op/\-\-output_path                 | pathway for output files to be written (default: same as -f input)                                                                           |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
| -rih, \-\-report_inparalog_handling | create a summary file of which inparalogs where kept compared to trimmed                                                                     |
+-------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------+
*For genome-scale analyses, we recommend changing the -o/\-\-occupancy parameter to be the same for all large gene families so that the minimum SNAP-OG occupancy is the same
for all SNAP-OGs.
