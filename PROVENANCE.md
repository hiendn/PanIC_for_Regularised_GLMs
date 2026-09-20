# Provenance

## Frozen implementation lineage

The confirmatory implementation inherits the scenarios, candidate paths,
cross-signed pilot targets, splits, grids, numerical controls, and
largest-multiplier one-standard-error rule from a preceding locked release.
The three inherited sources are vendored in `reference_implementation/` and
retain their production SHA-256 hashes:

```text
537d3cd266dd880453d194de9c03a7cbbbcdbaa7352f0b197886818d0dd97a5e  reference_implementation/config.R
9128cd0112cd3c679fef5d5601c22a4232293e2c38a9dae7e7caa6d2f3762c68  reference_implementation/path_functions.R
5df27ca79d5da1e7440b6c7b23874bdbdf62772fd02680fe9822136be8ffe568  reference_implementation/revised_cf_functions.R
```

The only primary-method amendment was the asymmetric undershoot weight

```text
omega(n_validation) = sqrt(log(log(n_validation + exp(exp(1))))).
```

`PanIC-CF-original` evaluates the preceding log-log weight on the identical
paths and targets. `CV-min` is the unchanged ordinary five-fold minimum-loss
comparator; `CV-1SE` is calculated from the same fold-loss curves and adds no
fit.

## Development evidence and untouched confirmation

The square-root log-log choice was motivated by a separate 100-replication
development experiment using master seed `2046091802`. Neither that output nor
the earlier 500-replication result was pooled with the confirmation.

The method, primary estimands, 0.1% prediction margin, sensitivity margins, and
decision rules were recorded before production outcomes were inspected. Early
two-replication packaging smoke runs consumed the then-candidate seeds
`2056091802` and `2076091802`; no statistical setting was changed in response.
Those seeds were retired so that the final production seeds were literally
untouched. The retained main and grid studies use `2066091802` and
`2086091802`; dedicated smoke mode uses `2106091802` and `2126091802`.

The deterministic audit enumerates the training, five calibration-split,
CV-fold, and test streams. It verifies exact pairwise non-overlap among all
four retained production/smoke families and conservative supersets of the
earlier or retired `2026091801`, `2026091802`, `2036091802`, `2046091802`,
`2056091802`, and `2076091802` families.

The original production method-lock file had SHA-256
`a7dcbff623c335d80a2a5dc710282d8fbd83c1f0b4ce3613e7f470fcd506db26`.
The repository manifest has different path-level hashes only because the
candidate filenames and relative source locations were normalized for a
self-contained public checkout. The statistical expressions and simulation
logic are unchanged. `METHOD_LOCK_SHA256.txt` records and verifies the current
repository paths.

## Production execution

The retained main run used seven settings and 1,000 replications per setting;
the grid study used 500 common-random-number replications at 61, 121, and 241
radii. The production environment recorded R 4.4.0, `glmnet` 4.1-10, and
`Matrix` 1.7-0 on `aarch64-apple-darwin20`.

The main run completed 7,000 data sets, 35,000 primary method rows, and 70,000
calibration-direction rows. It recorded zero method, full-path, calibration,
and CV failures; zero default uses; zero solver warnings; zero multiplier
endpoint selections; and zero PanIC-CF/CV selected-radius endpoint selections.
There was one projected target among 70,000 calibration rows. The maximum
radius interpolation error was `4.0274e-08`, below the locked `5e-08`
tolerance.

An independent verifier reconstructed the production summaries from the raw
RDS objects and matched the retained summaries to less than `5e-16`. All 21
deterministic checks, the full main result verifier, and the grid verifier
passed. The full result directory's file names, byte sizes, and SHA-256 hashes
at completion of that run are retained in
`results/production_file_manifest.csv` even though the bulky replication-level
files are not committed. This is intentionally an immutable historical
snapshot and therefore retains the original `manuscript_generated` paths. The
later self-contained release archive supplies a separate manifest for its
actual payload.

## Repository packaging

Repository packaging made only auditable, non-statistical changes:

- the frozen preceding sources were copied into `reference_implementation/`;
- scripts received descriptive public filenames and relative references;
- the default output directory was normalized from `results_v2` to `results`;
- compact production summaries, the exact seed ledger, and manuscript-facing
  tables were retained while regenerable raw rows were ignored;
- the runtime figure was re-rendered from the archived timings with the neutral
  labels `PanIC-CF`, `Shared 5-fold CV`, and `BIC-like`;
- reader-facing tables consistently label nonlinear active-count rows
  `BIC-like` while retaining the explicit exploratory-scope qualification.

No production outcome was recomputed or selected during packaging.

## Self-contained-table and denominator-audit amendment

The subsequent self-contained-table release made a bounded presentation and
reproducibility amendment. `manuscript/main.tex` now embeds all seven final
table environments directly between stable marker comments, so it has no
external table-file dependency. The files under
`results/table_*.tex` remain deterministic, checksum-controlled mirrors for
regeneration and audit. `Manuscript_Table_Rendering.R` and
`Render_Manuscript_Tables.R` reconstruct the seven mirrors from the retained
compact summaries, while `Manuscript_Table_Tools.R` refreshes or verifies the
embedded blocks byte-for-byte and rejects any remaining external table
reference.

The same release makes the domain of the paired relative test-deviance
estimand explicit. Every `CV-min` denominator must be finite and strictly
positive; analysis and independent verification fail on any violation. No
replication is discarded, and no floor or replacement is applied. The added
`relative_deviance_denominator_audit.csv` records the realised denominator
range and finite-sample relative-contrast moments by scenario, and the
full-result verifier reconstructs it from the archived paired rows.

These changes do not alter the statistical procedure, estimands, margins,
seeds, simulation rows, retained summaries, or decision rule. Rendering the
tables and reconstructing the denominator audit from the retained results
leave every reported numerical result unchanged.
