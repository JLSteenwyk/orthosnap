# Repository Guidelines

## Project Structure & Module Organization
Core source lives in `orthosnap/`, with the CLI entry point in `orthosnap/orthosnap.py` and shared utilities split across helper modules (for example `helper.py`, `writer.py`). Tests mirror runtime paths under `tests/` (`unit/` for isolated logic, `integration/` for CLI flows) and reuse fixtures in `tests/samples/`. Documentation is authored with Sphinx in `docs/`; build artifacts land in `docs/_build/` during CI. Reusable data for demonstrations ships in `samples/`.

## Build, Test, and Development Commands
Create an isolated environment (`python -m venv .venv && source .venv/bin/activate`) before installing. Run `make install` to install the package locally so the `orthosnap` console script is available. Use `make test.fast` for the default CI-equivalent quick check, or `make test.unit` and `make test.integration` when focusing on a specific layer. `make test.coverage` produces `unit.coverage.xml` and `integration.coverage.xml` for upload to Codecov; clean any generated FASTA outputs before committing.

## Coding Style & Naming Conventions
Target Python 3.9–3.13 and follow PEP 8: four-space indentation, snake_case functions, PascalCase classes, and descriptive module-level constants. Keep functions small, prefer explicit variable names (e.g., `taxa_counts` over `tc`), and include short docstrings for public helpers. Align logging and user-visible messages with existing phrasing in `orthosnap/helper.py`. When adding CLI arguments, route parsing through `parser.py` to remain consistent.

## Testing Guidelines
Pytest drives all suites. Mark unit-only tests with the default markers and integration cases with `@pytest.mark.integration`. Place new unit tests beside the feature under `tests/unit/` and use fixtures from `tests/samples/` where possible. Run `python -m pytest -k <pattern>` to iterate during development, and ensure `make test` succeeds before opening a pull request. Aim to keep coverage steady; add regression tests for every CLI flag or branch condition.

## Commit & Pull Request Guidelines
Write concise, present-tense commit messages (e.g., `add delimiter argument`) and group related changes. Reference issues in the body using `Fixes #<id>` when relevant. For pull requests, summarise behaviour changes, call out new dependencies, and attach command outputs or screenshots when UX changes. Confirm CI passes, ensure docs stay accurate (update Sphinx sources when CLI flags change), and request at least one review before merging.

## Environment & Documentation Tips
Use `requirements.txt` for runtime deps (supports CPython 3.9 through 3.13) and install doc tooling via `pipenv` inside `docs/` when updating the Sphinx site (`cd docs && pipenv run make html`). Avoid checking in files under `docs/_build/` or temporary trees left by integration tests.

## Performance Improvements In Progress
- Replace whole-tree `copy.deepcopy` calls in `orthosnap/orthosnap.py` with subtree-level cloning so we only touch the relevant clade per iteration.
- Precompute a parent lookup for tip nodes to let `get_subtree_tips` reuse cached structure instead of copying the tree for every duplicate gene check.
- Track `assigned_tips` as a set throughout execution to skip repeated list-to-set conversions while filtering already-handled sequences.
