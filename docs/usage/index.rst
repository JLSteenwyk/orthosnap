Usage
=====

This section covers practical OrthoSNAP usage. For a full worked example, see the tutorial.

OrthoSNAP takes:

- a gene tree in Newick format
- the FASTA file used to infer that tree

It outputs one FASTA file per inferred SNAP-OG (single-copy orthologous subgroup). Optionally, it can also write a Newick tree per SNAP-OG, an inparalog handling report, and one color-coded subgroup plot for the full input tree.

Basic usage
-----------

For most cases, only `-f/--fasta` and `-t/--tree` are required:

.. code-block:: shell

   $ orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre

Input requirements
------------------

- FASTA headers and tree tip labels must match.
- Taxon and sequence IDs must be separated by the same delimiter in both files.
- Default delimiter is `|` (for example, `species_A|gene_001`).

Accounting for tree uncertainty
-------------------------------

OrthoSNAP can collapse low-support bipartitions before pruning inparalogs.

- Default support threshold is `80`.
- Use `-s/--support` to change it.

.. code-block:: shell

   $ orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre -s 70

Choosing which inparalog to keep
--------------------------------

Use `-ip/--inparalog_to_keep` to select how species-specific inparalogs are resolved.

Supported values:

- `shortest_seq_len`
- `median_seq_len`
- `longest_seq_len` (default)
- `shortest_branch_len`
- `median_branch_len`
- `longest_branch_len`

Examples:

.. code-block:: shell

   $ orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre -ip shortest_branch_len

.. code-block:: shell

   $ orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre -ip median_seq_len

Inparalog handling report
-------------------------

Use `-rih/--report_inparalog_handling` to write a tab-delimited report named
`<input_fasta>.inparalog_report.txt`.

Columns are:

- SNAP-OG identifier
- kept inparalog
- trimmed inparalog(s), separated by `;`

.. code-block:: shell

   $ orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre -rih

Specifying the delimiter
------------------------

If your headers do not use `|`, specify the delimiter with `-d/--delimiter`.

.. code-block:: shell

   $ orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre -d -

Plotting SNAP-OG assignments
----------------------------

Use `-ps/--plot_snap_ogs` to create one figure of the full tree with distinct colors for each inferred SNAP-OG.
Default plot format is PNG; choose PDF or SVG with `-pf/--plot_format`.

.. code-block:: shell

   $ orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre -ps

.. code-block:: shell

   $ orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre -ps -pf svg

Example output (`png`):

.. image:: ../_static/img/orthosnap_subgroups_example.png
   :width: 100%
   :align: center
   :alt: Example OrthoSNAP subgroup plot with color-coded SNAP-OG assignments on a phylogeny.

All options
-----------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Option
     - Meaning
   * - ``-h/--help``
     - Print help message.
   * - ``-v/--version``
     - Print software version.
   * - ``-f/--fasta``
     - Input FASTA file.
   * - ``-t/--tree``
     - Input tree file in Newick format.
   * - ``-s/--support``
     - Collapse threshold for branch support (default: 80).
   * - ``-o/--occupancy``
     - Minimum represented taxa for subgroup candidates (default: rounded half of taxa in input FASTA).
   * - ``-r/--rooted``
     - Treat input tree as rooted; otherwise midpoint-root it (default: false).
   * - ``-d/--delimiter``
     - Delimiter between taxon and sequence IDs (default: ``|``).
   * - ``-st/--snap_trees``
     - Also write SNAP-OG trees in Newick format (default: false).
   * - ``-ip/--inparalog_to_keep``
     - Rule for keeping one inparalog among species-specific duplicates (default: ``longest_seq_len``).
   * - ``-rih/--report_inparalog_handling``
     - Write tab-delimited inparalog handling report (default: false).
   * - ``-op/--output_path``
     - Output directory (default: directory containing input FASTA).
   * - ``-ps/--plot_snap_ogs``
     - Write one color-coded full-tree plot with subgroup labels (default: false).
   * - ``-pf/--plot_format``
     - Output format for subgroup plot (``png`` default, or ``pdf``/``svg``).

For genome-scale analyses, consider using the same `-o/--occupancy` value across all gene families to keep SNAP-OG occupancy thresholds consistent.
