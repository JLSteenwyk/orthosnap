import sys
import textwrap

from argparse import (
    ArgumentParser,
    RawTextHelpFormatter,
    SUPPRESS,
    RawDescriptionHelpFormatter,
)

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
        <tree>                                      asdf

        <fasta>                                     asdf                 
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
        -s, --support                                  support
        <support>
        
        -------------------------------------
        | Detailed explanation of arguments | 
        -------------------------------------
        -s, --support <support>
            support to collapse tree at
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
        "-h",
        "--help",
        action="help",
        help=SUPPRESS,
    )

    optional.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"clipkit {__version__}",
        help=SUPPRESS,
    )

    return parser
