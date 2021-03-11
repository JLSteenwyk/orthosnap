import textwrap
import time


def write_user_args(tree: str, fasta: str, support: float, occupancy: float):
    """
    Function to print user arguments to stdout
    """
    print(
        textwrap.dedent(
            f"""\
    -------------
    | Arguments |
    -------------
    Input phylogeny: {tree}
    Input fasta: {fasta}
    Support threshold: {support}
    Taxon occupancy threshold: {occupancy}
    """
        )
    )


def write_output_stats(fasta, subgroup_counter, start_time):
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
