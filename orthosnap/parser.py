import textwrap

from argparse import (
    ArgumentParser,
    SUPPRESS,
    RawDescriptionHelpFormatter,
)

from .helper import InparalogToKeep

from .version import __version__


def create_parser():
    parser = ArgumentParser(
        add_help=False,
        formatter_class=RawDescriptionHelpFormatter,
        usage=SUPPRESS,
        description=textwrap.dedent(
            f"""\
        OrthoSNAP identifies single-copy orthologous subgroups (SNAP-OGs)
        from larger multi-copy gene families using a gene tree and FASTA file.

        Version: {__version__}
        Citation: Steenwyk et al. (2022), PLOS Biology.
        DOI: 10.1371/journal.pbio.3001827
        https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001827

        Usage: orthosnap -f <fasta_file> -t <newick_tree_file> [options]
        Example:
          orthosnap -f OG0000010.fa -t OG0000010.treefile -s 80 -o 5
        """
        ),
    )

    ## required arguments
    required = parser.add_argument_group(
        "required arguments",
        description=textwrap.dedent(
            """\
        -t, --tree <newick tree file>
            Input phylogeny in Newick format.

        -f, --fasta <fasta file>
            Input sequences in FASTA format.
        """
        ),
    )

    required.add_argument(
        "-f",
        "--fasta",
        type=str,
        required=False,
        help=SUPPRESS,
        metavar="fasta",
    )

    required.add_argument(
        "-t",
        "--tree",
        type=str,
        required=False,
        help=SUPPRESS,
        metavar="tree",
    )

    ## optional arguments
    optional = parser.add_argument_group(
        "optional arguments",
        description=textwrap.dedent(
            """\
        -s, --support  <support>
            Collapse bipartitions with support values below this threshold.
            Default: 80

        -o, --occupancy  <occupancy>
            Minimum number of taxa required for a subgroup candidate.
            Default: rounded half of total taxa in the input FASTA

        -r, --rooted
            Treat the input tree as already rooted.
            Default: false (tree is midpoint-rooted)

        -d, --delimiter
            Delimiter between taxon and sequence name in headers.
            Default: "|"

        -st, --snap_trees
            Also write Newick trees for inferred SNAP-OGs.
            Default: false

        -ip, --inparalog_to_keep <shortest_seq_len,
                                  median_seq_len,
                                  longest_seq_len,
                                  shortest_branch_len,
                                  median_branch_len,
                                  longest_branch_len>
            Rule used to keep one sequence among species-specific inparalogs.
            Default: longest_seq_len

        -rih, --report_inparalog_handling
            Write a report of kept vs trimmed inparalogs.

        -op, --output_path
            Directory for output files.
            Default: same directory as input FASTA

        -ps, --plot_snap_ogs
            Plot full tree with color-coded SNAP-OG assignments.
            Default: false

        -pf, --plot_format <png|pdf|svg>
            Output format for SNAP-OG plot.
            Default: png

        Notes
        -----
        -t, --tree <newick tree file>
            Input tree in Newick format.

        -f, --fasta <fasta file>
            Input sequence file in FASTA format. Header naming must match the tree.
            Example header format: species_A|gene_001
        
        -s, --support  <support>
            Typical values:
            80 for ultrafast bootstrap support; 70 for standard bootstrap support.

        -o, --occupancy  <occupancy>
            Occupancy can be any positive float.
            Subgroups are considered when represented taxa >= occupancy.

        -r, --rooted
            If omitted, OrthoSNAP midpoint-roots the tree.

        -d, --delimiter
            Must appear in all relevant FASTA/tree labels.
            Example: >species_A|gene0

        -st, --snap_trees
            Writes additional .tre files for each output SNAP-OG.

        -ip, --inparalog_to_keep <shortest_seq_len,
                                  median_seq_len,
                                  longest_seq_len,
                                  shortest_branch_len,
                                  median_branch_len,
                                  longest_branch_len>                                
            Sequence-based options use ungapped sequence length.
            Branch-based options use tip-to-root distance in the current subtree.

        -rih, --report_inparalog_handling
            Report columns:
            1) SNAP-OG identifier
            2) kept inparalog
            3) trimmed inparalog(s), separated by ';'

        -op, --output_path <str>
            Parent directory for output files.

        -ps, --plot_snap_ogs
            Writes one color-coded tree figure showing inferred SNAP-OGs.

        -pf, --plot_format <png|pdf|svg>
            File format for the SNAP-OG assignment plot.
        """
        ),
    )

    optional.add_argument(
        "-s",
        "--support",
        type=float,
        required=False,
        help=SUPPRESS,
        metavar="seed",
    )

    optional.add_argument(
        "-o",
        "--occupancy",
        type=float,
        required=False,
        help=SUPPRESS,
        metavar="seed",
    )

    optional.add_argument(
        "--occupancy-count",
        type=int,
        required=False,
        help=SUPPRESS,
        metavar="count",
    )

    optional.add_argument(
        "--occupancy-fraction",
        type=float,
        required=False,
        help=SUPPRESS,
        metavar="fraction",
    )

    optional.add_argument(
        "-r",
        "--rooted",
        action="store_true",
        required=False,
        help=SUPPRESS,
    )

    optional.add_argument(
        "-d",
        "--delimiter",
        type=str,
        required=False,
        help=SUPPRESS,
    )

    optional.add_argument(
        "-st",
        "--snap_trees",
        action="store_true",
        required=False,
        help=SUPPRESS,
    )

    inparalog_to_keep_choices = [val.value for val in InparalogToKeep]
    optional.add_argument(
        "-ip",
        "--inparalog_to_keep",
        help=SUPPRESS,
        required=False,
        nargs="?",
        choices=inparalog_to_keep_choices,
    )

    optional.add_argument(
        "-rih",
        "--report_inparalog_handling",
        action="store_true",
        required=False,
        help=SUPPRESS,
    )

    optional.add_argument(
        "-op",
        "--output_path",
        help=SUPPRESS,
        required=False,
    )

    optional.add_argument(
        "--manifest",
        type=str,
        required=False,
        help=SUPPRESS,
    )

    optional.add_argument(
        "--validate-only",
        action="store_true",
        required=False,
        help=SUPPRESS,
    )

    optional.add_argument(
        "--resume",
        action="store_true",
        required=False,
        help=SUPPRESS,
    )

    optional.add_argument(
        "--structured-output",
        action="store_true",
        required=False,
        help=SUPPRESS,
    )

    optional.add_argument(
        "--bootstrap-trees",
        type=str,
        required=False,
        help=SUPPRESS,
    )

    optional.add_argument(
        "--consensus-min-frequency",
        type=float,
        required=False,
        help=SUPPRESS,
        metavar="fraction",
    )

    optional.add_argument(
        "-ps",
        "--plot_snap_ogs",
        action="store_true",
        required=False,
        help=SUPPRESS,
    )

    optional.add_argument(
        "-pf",
        "--plot_format",
        type=str,
        choices=["png", "pdf", "svg"],
        required=False,
        help=SUPPRESS,
    )

    optional.add_argument(
        "-h",
        "--help",
        action="help",
        help=SUPPRESS,
    )

    optional.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"orthosnap {__version__}",
        help=SUPPRESS,
    )

    return parser
