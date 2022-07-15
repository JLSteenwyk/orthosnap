import textwrap
import time

from orthosnap.helper import InparalogToKeep


def write_user_args(
    tree: str,
    fasta: str,
    support: float,
    occupancy: float,
    rooted: bool,
    snap_trees: bool,
    inparalog_to_keep: InparalogToKeep,
):
    """
    Function to print user arguments to stdout
    """
    print(
        textwrap.dedent(
            f"""\
    -------------
    | Arguments |
    -------------
    Input phylogeny: {tree} (rooted, {rooted})
    Input fasta: {fasta}
    Support threshold: {support}
    Taxon occupancy threshold: {occupancy}
    Output newick of SNAP-OGs: {snap_trees}
    Inparalog to keep: {inparalog_to_keep.value}
    """
        )
    )


def write_output_stats(fasta, subgroup_counter, start_time, snap_trees):
    """
    Function to print out output statistics
    """
    if subgroup_counter > 0:
        print(
            textwrap.dedent(
                f"""\

            ---------------------
            | Output Statistics |
            ---------------------
            Subgroups of single-copy orthologous genes identified: {subgroup_counter}
            Output files:"""
            )
        )
        if snap_trees:
            for i in range(subgroup_counter):
                print(f"\t{fasta}.orthosnap.{i}.fa\n\t{fasta}.orthosnap.{i}.tre")
        else:
            for i in range(subgroup_counter):
                print(f"\t{fasta}.orthosnap.{i}.fa")
            print(
                textwrap.dedent(
                    f"""\
                Execution time: {round(time.time() - start_time, 3)}s
        """
                )
            )
    else:
        print(
            textwrap.dedent(
                f"""\

            ---------------------
            | Output Statistics |
            ---------------------
            Single-copy orthologous genes identified: {subgroup_counter}
            Execution time: {round(time.time() - start_time, 3)}s
        """
            )
        )
