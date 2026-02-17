import copy
import re

from Bio import Phylo


def _iter_clade_descendants(clade):
    stack = [clade]
    while stack:
        current = stack.pop()
        yield current
        stack.extend(current.clades)


def _to_hex_color(rgba):
    return "#{:02x}{:02x}{:02x}".format(
        int(rgba[0] * 255),
        int(rgba[1] * 255),
        int(rgba[2] * 255),
    )


def plot_snap_ogs(
    tree,
    subgroup_records: list,
    fasta: str,
    output_path: str,
    plot_format: str,
):
    """
    Plot full phylogeny and color SNAP-OG subgroup clades and tips.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError as exc:
        raise ImportError(
            "matplotlib is required for --plot_snap_ogs. Install it in your environment."
        ) from exc

    if not subgroup_records:
        return None

    tree_for_plot = copy.deepcopy(tree)

    cmap = plt.get_cmap("tab20")
    subgroup_colors = dict()
    tip_colors = dict()

    ordered_records = sorted(subgroup_records, key=lambda rec: rec["subgroup_id"])
    for idx, record in enumerate(ordered_records):
        color = _to_hex_color(cmap(idx % cmap.N))
        subgroup_colors[record["subgroup_id"]] = color
        tips = record["tips"]
        for tip in tips:
            tip_colors[tip] = color
        try:
            subgroup_clade = tree_for_plot.common_ancestor(tips)
        except Exception:
            continue
        for descendant in _iter_clade_descendants(subgroup_clade):
            descendant.color = color

    terminal_count = len(tree_for_plot.get_terminals())
    fig_height = max(8, min(0.22 * terminal_count, 28))
    fig, ax = plt.subplots(figsize=(14, fig_height))

    Phylo.draw(
        tree_for_plot,
        axes=ax,
        do_show=False,
        label_func=lambda clade: clade.name if clade.is_terminal() else None,
    )

    for text in ax.texts:
        label = text.get_text().strip()
        if label in tip_colors:
            text.set_color(tip_colors[label])
            text.set_fontweight("bold")

    legend_handles = []
    for subgroup_id in sorted(subgroup_colors):
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color=subgroup_colors[subgroup_id],
                lw=3,
                label=f"SNAP-OG {subgroup_id}",
            )
        )
    ax.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=max(1, min(4, len(legend_handles))),
        frameon=False,
        title="Subgroups",
    )

    ax.set_ylabel("")
    ax.tick_params(axis="y", which="both", left=False, labelleft=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    fig.tight_layout(rect=[0.0, 0.08, 1.0, 1.0])

    fasta_path_stripped = re.sub("^.*/", "", fasta)
    output_file = (
        f"{output_path}{fasta_path_stripped}.orthosnap.subgroups.{plot_format}"
    )

    save_kwargs = {"format": plot_format, "bbox_inches": "tight"}
    if plot_format == "png":
        save_kwargs["dpi"] = 300
    fig.savefig(output_file, **save_kwargs)
    plt.close(fig)

    return output_file
