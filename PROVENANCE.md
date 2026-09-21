# Provenance

## Frozen numerical lineage

The project retains the preceding generic numerical engine under
`reference_implementation/`. Its three files keep their historical SHA-256
hashes:

```text
537d3cd266dd880453d194de9c03a7cbbbcdbaa7352f0b197886818d0dd97a5e  reference_implementation/config.R
9128cd0112cd3c679fef5d5601c22a4232293e2c38a9dae7e7caa6d2f3762c68  reference_implementation/path_functions.R
5df27ca79d5da1e7440b6c7b23874bdbdf62772fd02680fe9822136be8ffe568  reference_implementation/revised_cf_functions.R
```

These files preserve the inherited scenarios, path construction, split logic,
grids, and numerical controls. Current top-level files may source their generic
primitives, but the frozen directory is not the active assessed implementation:
it does not determine the current method set, exact support rule, seeds, or
output schema.

## Current locked revision

The active revision assesses PanIC-CF, ordinary minimum-loss five-fold `CV`,
and the family-appropriate `BIC-like` comparator. It also replaces the inherited
coefficient threshold with the literal fitted-nonzero rule for both support
evaluation and BIC-like active counts.

The three-method, exact-zero procedure is unchanged in the present revision.
The primary design now crosses five data-generating designs with both
`n = 500` and `n = 1000`, yielding ten settings. Earlier exploratory and
historical results are not pooled with the regenerated results.

The square-root log-log calibration weight remains

```text
omega(n_validation) = sqrt(log(log(n_validation + exp(exp(1))))).
```

Each setting retains 1,000 replications, a 2,000-observation test sample, the
same radius and multiplier grids, and five bidirectional half-splits. The
pooled PanIC-CF/CV calculations and the `0.001` reference margin remain in an
internal diagnostic file, but no pooled pass condition is a release gate.

The auxiliary design uses boundary sample sizes `n = 500` and `n = 1000`.
Its subsampling illustration uses a full sample of `n = 1000` and subsample
sizes `b = 50, 100, 200`.

## Seed separation

The active families are:

- main production: `2146092101`;
- grid production: `2138092001`;
- runtime benchmark: `2096092101`;
- main smoke: `2140092001`; and
- grid smoke: `2142092001`.

The configuration treats the following earlier families as prior or retired:

```text
2026091801  2026091802  2036091802  2046091802  2056091802
2066091802  2076091802  2086091802  2106091802  2126091802
2136092001
```

This list includes the retired main seed `2136092001`, earlier main, grid, and
dedicated smoke families, and development candidates. The deterministic audit
enumerates training, five calibration-split, CV-fold, and test streams and
checks the new production, smoke, and runtime families for pairwise non-overlap,
non-overlap with conservative supersets of every prior family, and validity in
R's positive integer seed range.

Smoke mode is explicit, uses only its dedicated seeds, rejects the committed
`results/` directory, and cannot consume a production seed. Reduced-replication
runs likewise cannot masquerade as a locked production run.

## Regeneration and verification

The main regeneration produces 10,000 data sets and 30,000 primary method rows:
one row for each of the three assessed methods in every replication. The grid
regeneration uses 500 common-random-number replications at 61, 121, and 241
radii. The runtime benchmark also uses the same three reader-facing method
labels.

The main and grid runners record their configurations, environments, seed
roles, explicit ledgers, raw rows, diagnostics, and compact summaries.
Independent scripts reconstruct the main and grid summaries and the internal
pooled diagnostic from the replication-level outputs.
`confirmatory_decision.csv` is retained only for that descriptive audit; its
legacy-named pass fields do not determine release status. Release checksums and
manifests are refreshed after the new runs, so this document does not carry
forward numerical claims from a superseded seed family.

The intended reproduction environment is R 4.4.0 with `glmnet` 4.1-10 and
`Matrix` 1.7-0 on `aarch64-apple-darwin20`. The recorded environment files,
rather than this prose alone, are authoritative for an executed run.

## Repository packaging

Packaging keeps the statistical and presentation layers distinct:

- current top-level scripts define and verify the active three-method procedure;
- `reference_implementation/` remains frozen for lineage and numerical audit;
- replication-level files are regenerable, while compact summaries, table
  mirrors, configurations, ledgers, and verification records form the release
  record; and
- manuscript labels are consistently `PanIC-CF`, `CV`, and `BIC-like`.

The manuscript embeds six complete table environments between stable marker
comments. The deterministic, checksum-controlled mirrors are the primary
support, primary performance, calibration, grid-sensitivity, boundary, and
subsampling tables; no confirmatory-decision table is generated. The primary
performance table omits the redundant attained-radius column, the calibration
table leaves sparse operational counts to the diagnostics and prose, and the
grid table omits runtime. Those layout choices do not change an estimand or
fitted result.

The relative-deviance audit enforces the endpoint's domain: every paired CV
denominator must be finite and strictly positive. A violation is a hard failure;
no observation is deleted and no denominator is replaced.
