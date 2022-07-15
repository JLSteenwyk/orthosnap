import sys
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
                   _   _                                 
                  | | | |                                
         ___  _ __| |_| |__   ___  ___ _ __   __ _ _ __  
        / _ \| '__| __| '_ \ / _ \/ __| '_ \ / _` | '_ \ 
       | (_) | |  | |_| | | | (_) \__ \ | | | (_| | |_) |
        \___/|_|   \__|_| |_|\___/|___/_| |_|\__,_| .__/ 
                                                  | |    
                                                  |_|     

        Version: {__version__}
        Citation: Steenwyk et al.
        Usage: orthosnap <input> [optional arguments]
        """
        ),
    )

    # if no arguments are given, print help and exit
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()

    ## required arguments
    required = parser.add_argument_group(
        "required arguments",
        description=textwrap.dedent(
            """\
        -t, --tree <newick tree file>
            input tree file in newick format

        -f, --fasta <fasta file>
            input sequence file in fasta format                 
        """
        ),
    )

    required.add_argument(
        "-f",
        "--fasta",
        type=str,
        required=True,
        help=SUPPRESS,
        metavar="fasta",
    )

    required.add_argument(
        "-t",
        "--tree",
        type=str,
        required=True,
        help=SUPPRESS,
        metavar="tree",
    )

    ## optional arguments
    optional = parser.add_argument_group(
        "optional arguments",
        description=textwrap.dedent(
            """\
        -s, --support  <support>
            support threshold for bipartition collapsing
            default: 80

        -o, --occupancy  <occupancy>
            occupancy threshold for minimum number of tips in orthologous subgroup
            default: 50 percent of total number of taxa  
       
        -r, --rooted
            boolean argument for whether the input phylogeny is already rooted
            default: false

        -st, --snap_trees
            boolean argument for whether trees of SNAP-OGs should be outputted
            default: false

        -ip, --inparalog_to_keep <shortest_seq_len,
                                  median_seq_len,
                                  longest_seq_len,
                                  shortest_branch_len,
                                  median_branch_len,
                                  longest_branch_len>
            determine which sequence to keep in the case of species-specific
            inparalogs using sequence- or tree-based options
            default: longest_seq_len

        
        -------------------------------------
        | Detailed explanation of arguments | 
        -------------------------------------
        -t, --tree <newick tree file>
            - input tree file in newick format
            - taxa name and gene name should be separate by a "|" symbol
              For example, "gene_a" from "species_a" should appear in the
              tree as "species_a|gene_a", and "gene_b" from "species_a"
              should appear in the tree as "species_a|gene_b", and so 
              on and so forth.

        -f, --fasta <fasta file>
            - input sequence file in fasta format 
            - taxa name and gene name should be formatted the same as in
              the tree file. Thus, "gene_a" from "species_a" should appear
              in the tree as "species_a|gene_a" and so on and so forth.
        
        -s, --support  <support>
            - support threshold for bipartition collapsing
            - all bipartitions will values less than the specified value 
              will be collapsed. For example, if the support threshold
              value is 80, all bipartitions with 79 or less support
              will be collapsed.
            - default value is 80 and is set for ultrafast bootstrap
              approximations. If bipartitions support was evaluated 
              using standard bootstrap, a common threshold to use is 70.

        -o, --occupancy  <occupancy>
            - occupancy threshold for minimum number of tips in orthologous
              subgroup. 
            - default value is 50 percent of the total number of taxa
            - values are rounded to the nearest integer. For example,
              if there are 15 taxa, the occupancy threshold will be 8.

        -r, --rooted
            - boolean argument for whether the input phylogeny is already rooted
            - if used, the input phylogeny is assumed to be rooted; if not,
              the tree will be midpoint rooted

        -st, --snap_trees
            - boolean argument for whether newick files of SNAP-OGs should also
              be outputted
            - if used, newick files of SNAP-OGs will be outputteds

        -ip, --inparalog_to_keep <shortest_seq_len,
                                  median_seq_len,
                                  longest_seq_len,
                                  shortest_branch_len,
                                  median_branch_len,
                                  longest_branch_len>                                
            - specify how to determine which species-specific inparalog should be kept
            - the species-specific inparalog can be kept based on sequence length
              (shortest/median/longest_seq_len) or branch length based on tip-to-root
              distances (shortest/median/longest_branch_len)
            - by default, the longest sequence is kept following the standard approach
              in transcriptomics
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
        "-r",
        "--rooted",
        action="store_true",
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
