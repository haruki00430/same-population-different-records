# Analysis audit

Generated from executed results at 2026-09-19T23:14:37.139126+00:00.

## Decision: GO — strong under the predefined illustrative criterion

The prespecified 5-percentage-point ascertainment change (0.80 to 0.85) produces an analytical relative measure change of 6.250% with clinical event risk fixed. Monte Carlo mean relative change is 6.354%. The 10% safety improvement masking root is 8.889 percentage points. These meet the operational <=10 pp criterion. This does not establish that the perturbation is empirically common or clinically plausible in any particular setting.

## Validation and completeness

- 1085 A/C cells checked against analytical expectations; maximum absolute Monte Carlo z discrepancy = 3.564; required <6.
- Main A/C cells: 409; main risk cells: 168; sensitivity A/C cells: 676; sensitivity risk cells: 88; SPC rule rows: 80 (two rules per scenario).
- Every listed simulation cell has 10,000 replicates. All grid cells, including unchanged records and weak denominator effects, are retained.
- Protocol/config hashes and execution log document this run. No empirical patient data were used. Development model fit gradient infinity norm = 7.93e-11.
- B1 baseline O/E is not assumed calibrated. B2 is fixed after development; development-sample uncertainty is not propagated.
- A/C identity and analytical tests passed before full execution. Full-run verification is saved separately in audit/test_results.txt.

## Results that constrain interpretation

- A 20% clinical event reduction is erased only at perfect ascertainment (1.00) from baseline .80; population-direction reversal is unattainable through sensitivity alone in this main model (finite-sample direction still varies).
- Non-differential denominator sensitivity cancels exactly when denominator specificity is perfect. Nonzero FPR breaks that cancellation, but effects may remain small. This negative control is retained.
- Tipping points for risk models outside c_B=0.50–1.00 are labelled as outside the main grid, not presented as observed grid reversals.
- Signal probability is compared with no-change SPC controls; an individual signal cannot be causally attributed to recording on the basis of the chart alone.
- Paired probabilities depend on the recording-error coupling. Independent-cohort and independent-recording analyses are reported.
- Empirical replicate intervals, Monte Carlo uncertainty, nominal sampling CIs and signal/noise reliability have distinct meanings.

## Editor's two-sentence question

Previous work establishes that recording practices can influence quality measures. This study provides an explicit counterfactual stress-test that holds clinical reality constant, maps changes in record generation to changes in quality inference, and identifies the recording changes sufficient to erase or reverse improvement across measure architectures.

## Submission checks

Verify current BMJ Quality & Safety article-type availability, word/figure limits and declarations manually; the official author page was blocked during this session. Add authors, affiliations, funding, contributions, conflicts and institutional ethics determination where applicable. No journal submission or protocol registration has occurred.
