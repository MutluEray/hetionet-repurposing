# Temporal validation plan (written before any results)

**Status: pre-registered, implemented in `src/13_temporal_validation.py`, not yet run on real data.** The pipeline is frozen at the commit tagged `frozen-v1`
(`git tag frozen-v1` before running anything below). Anything changed after seeing validation results is reported as
exploratory.

## Why

Hetionet v1.0 encodes knowledge available around 2015-16 (PharmacotherapyDB is dated March 2016; the network was
released Feb 2017; the paper appeared Sep 2017). Evidence that appeared afterwards can serve as a test set the
graph could not have seen. Rephetio itself validated on two external sets, so this follows precedent.

## What we test

Question: does the consensus score rank compounds that later entered clinical trials for the disease above what
simple baselines predict? It does **not** test efficacy; a trial entry only means someone found the hypothesis
worth testing, and negative trials count as positives.

* **Diseases:** Alzheimer's disease and Parkinson's disease (primary); pooled and separately.
* **Positives:** compounds present in Hetionet with at least one interventional drug trial in phase 2, 2/3 or 3 for the
  disease, start date on or after 2017-01-01 (ClinicalTrials.gov API v2), not already a Hetionet treats/palliates
  pair. Intervention names are matched to Hetionet compound names (exact, case-insensitive, salts stripped); every
  match is listed in the output for manual checking.
* **Ranking under test:** the un-triaged consensus score over all compounds (not the triaged shortlist, which used
  hand-made literature knowledge).
* **Baselines:** (1) number of genes the compound binds; (2) number of pre-2017 trials of the compound for the same
  disease (popularity of the hypothesis); (3) number of pre-2017 trials of the compound for any disease;
  (4) random. Comparator: Rephetio's published probabilities (het.io/repurpose), if available.
* **Primary metric:** AUROC of consensus minus AUROC of the best baseline, bootstrap 95% CI over positives.
  **Secondary:** hit@5%, average precision, and the same analysis restricted to compounds with no pre-2017 trial
  for that disease (novelty only).
* **Minimum evidence:** at least 20 pooled positives; otherwise the result is reported as inconclusive.

## Declared outcomes

| Outcome | What we would see | How we report it |
|---|---|---|
| Success | consensus beats every baseline, CI of the difference excludes 0, and the novelty-only analysis agrees | graph adds information beyond popularity |
| Partial | beats degree and random but not the pre-2017 trial-count baseline | graph mostly tracks existing research interest |
| Null | no better than degree or trial count | honest negative result; the benchmark gains came from popularity |
| Leakage | AUROC high only for compounds that already had pre-2017 trials | graph reflects literature, not foresight |
| Hub failure | good AUROC but top-20 dominated by promiscuous or endogenous compounds | passes globally, fails as a candidate generator |
| Direction failure | top candidates contraindicated or opposite-acting (needs a contraindication source, see below) | binding graph cannot support repurposing without direction |
| Instability | conclusions flip with w, permutation seed, or the disease | not robust |

I expect Partial or Null: trial entry is strongly driven by prior interest, and our earlier results show that part of
the benchmark signal is popularity. Any of these outcomes is reportable.

## Optional add-on: direction of effect

The largest known gap is that "binds" has no direction. A modern source with action types and contraindications
(for example PrimeKG, which has indication, contraindication and off-label drug-disease edges) could annotate the
candidates. Declared before use: it would be used only to measure and flag contraindicated or opposite-acting
candidates, not to tune the score, and its coverage of our candidates would be reported.

## Implementation notes and known deviations

* Trials come from the ClinicalTrials.gov API v2, cached with the retrieval date in `data/raw/ctgov_<slug>.json`.
* Name matching is exact after removing salts, doses and formulations, plus the alias list in
  `data/drug_synonyms.csv`. Fix matching only from the list of unmatched names (`ctgov_unmatched_*.csv`), never from
  outcomes, and report the first run and the final run.
* Only compounds that exist in Hetionet can be positives, so drugs introduced after about 2016 are invisible.
* The any-disease popularity baseline needs about 1,400 API calls (`--any-disease-counts`, roughly 30 minutes). If it
  is skipped, the report says so and that counts as a deviation.
* The Rephetio comparator (`--rephetio`) is best-effort: the script detects the compound, disease and probability columns
  and skips it if it cannot.
* Contraindication checking (the direction add-on) is not implemented in this script.
