.. _change_log:


Change log
==========

Major changes to OrthoSNAP are summarized here.

**1.6.0**
- Added batch manifest mode (`--manifest`) for running multiple orthogroups from one TSV/CSV.
- Added input preflight and dry-run validation mode (`--validate-only`).
- Added structured provenance outputs (`--structured-output`) with JSON and TSV summaries.
- Added explicit occupancy semantics (`--occupancy-count`, `--occupancy-fraction`).
- Added resumable execution (`--resume`) to skip already completed runs.
- Added bootstrap consensus subgrouping (`--bootstrap-trees`, `--consensus-min-frequency`, `--consensus-trees`).
- Improved CLI invocation handling by honoring `main(argv)` arguments.
- Relative to 1.5.0, this release emphasizes reproducibility and workflow orchestration rather than visualization.

**1.5.0**
- Added `-ps/--plot_snap_ogs` to generate one color-coded tree figure showing inferred SNAP-OG assignments.
- Added `-pf/--plot_format` to select `png` (default), `pdf`, or `svg` output for subgroup plots.
- Improved runtime for duplicate-handling hotspots via faster internal indexing and pruning updates.

**1.3.2**
Added argument for user-defined delimiter (-d, \-\-delimiter)

**1.2.0**
Added the -rih (\-\-report_inparalog_handling) function, which creates
a summary file of which inparalogs were kept compared to trimmed

**0.1.0**
Added -r/\-\-rooted, -st/\-\-snap_trees, and -ip/\-\-inparalog_to_keep functions

**1.0.0**
Improved species inparalog pruning
