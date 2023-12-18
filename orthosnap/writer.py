import re
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
    report_inparalog_handling: bool,
    output_path: str,
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
    Inparalog to keep: {inparalog_to_keep.value}
    Report inparalog handling: {report_inparalog_handling}
    Support threshold: {support}
    Taxon occupancy threshold: {occupancy}
    Output newick of SNAP-OGs: {snap_trees}
    Output directory: {output_path}
    """
        )
    )


def write_output_stats(
    fasta, subgroup_counter, start_time, snap_trees, output_path
):
    """
    Function to print out output statistics
    """

    fasta_path_stripped = re.sub("^.*/", "", fasta)
    output_file_name = (
        f"{output_path}/{fasta_path_stripped}.orthosnap.{subgroup_counter}.fa"
    )

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
                output_file_name = f"{output_path}{fasta_path_stripped}.orthosnap.{i}"
                print(f"\t{output_file_name}.fa\n\t{output_file_name}.tre")
        else:
            for i in range(subgroup_counter):
                output_file_name = f"{output_path}{fasta_path_stripped}.orthosnap.{i}"
                print(f"\t{output_file_name}.fa")
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
