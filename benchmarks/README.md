# OrthoSNAP performance benchmarks

`benchmark_performance.py` generates deterministic synthetic orthogroups in a
temporary directory. Each synthetic tree contains repeated taxon names across
otherwise single-copy subgroups, so it exercises candidate discovery,
duplicate handling, bootstrap aggregation, and consensus-tree pruning without
requiring a large fixture in the repository.

The benchmark supports balanced and pectinate (`unbalanced`) trees. Data
generation occurs outside measured regions. Runtime results use warm-up runs
followed by repeated measurements; memory results report peak Python
allocations measured by `tracemalloc`, not total process RSS. Every repeated run
must produce the same subgroup count and deterministic result signature.

## Reproduce

Run all benchmark modes with the default 640-tip pectinate tree:

```shell
python benchmarks/benchmark_performance.py --mode all
```

Use more repetitions for stable timing without the overhead of allocation
tracing:

```shell
python benchmarks/benchmark_performance.py \
  --mode core --shape unbalanced --groups 80 --taxa-per-group 8 \
  --warmups 2 --runs 7 --skip-memory
```

Measure larger scaling points:

```shell
python benchmarks/benchmark_performance.py \
  --mode core --shape unbalanced --groups 320 --taxa-per-group 8 \
  --warmups 1 --runs 3 --skip-memory
```

The `bootstrap` mode performs end-to-end parsing, extraction across five trees,
support aggregation, and consensus FASTA writing. The `consensus` mode isolates
consensus FASTA and induced-tree generation. The `startup` mode measures
`python -m orthosnap --help` in a fresh process.

## July 13, 2026 results

Environment: Apple arm64, macOS 26.4.1, Python 3.11.14, Biopython 1.85. The
baseline is commit `d12374c5`, which added this harness without changing runtime
code. The optimized implementation is commit `42a67e14`. Values are medians.

| Workload | Baseline | Optimized | Speedup |
| --- | ---: | ---: | ---: |
| Core, unbalanced, 80 groups × 8 taxa | 0.4731 s | 0.0103 s | 45.7× |
| Core, balanced, 80 groups × 8 taxa | 0.0641 s | 0.0039 s | 16.3× |
| Bootstrap, unbalanced, 20 groups × 8 taxa, 5 trees | 0.1381 s | 0.0140 s | 9.9× |
| Consensus trees, unbalanced, 40 groups × 8 taxa | 0.4906 s | 0.0262 s | 18.7× |
| CLI help startup | 0.2450 s | 0.2473 s | within noise |

| Workload | Baseline peak | Optimized peak | Reduction |
| --- | ---: | ---: | ---: |
| Core, unbalanced, 640 tips | 3.218 MiB | 0.465 MiB | 85.6% |
| Core, balanced, 640 tips | 1.314 MiB | 0.447 MiB | 66.0% |
| Bootstrap, 160 tips × 5 trees | 0.899 MiB | 0.579 MiB | 35.6% |
| Consensus trees, 320 tips | 1.281 MiB | 0.462 MiB | 63.9% |

Pectinate scaling improved as well. Doubling from 640 to 1,280 tips increased
baseline runtime from 0.4731 s to 2.3721 s (5.0×), while optimized runtime
increased from 0.0103 s to 0.0366 s (3.5×). The optimized 2,560-tip case
completed in 0.1322 s. Equivalent inputs retained identical benchmark result
signatures before and after optimization.

The measured gains come primarily from compact interval-based subtree metadata,
shared topology and parent indexes, logical inparalog pruning before any clone,
single-pass input reuse, batched FASTA serialization, and direct construction
of induced consensus trees.
