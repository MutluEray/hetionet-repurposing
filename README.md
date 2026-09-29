# Explainable knowledge-graph drug repurposing for Alzheimer's and Parkinson's disease (Hetionet)

Path-level, explainable repurposing hypotheses from a real 47k-node, 2.25M-edge biomedical knowledge graph
([Hetionet v1.0](https://github.com/hetio/hetionet), Himmelstein et al., eLife 2017) - with an honest evaluation
of what the graph can and cannot tell you.

**Not a discovery claim.** This is a methods and explainability project: it reimplements degree-weighted path
counting (the Project Rephetio approach), adds a permutation-null specificity score and gene-level counterfactuals,
and shows every candidate as the actual paths behind it. Full write-up: [`docs/writeup.md`](docs/writeup.md).

## Results at a glance

Pan-disease benchmark: rank every compound for each of 58 diseases with at least 3 known "treats" edges
(731 positives); macro AUROC, bootstrap 95% CI over diseases.

| Score | AUROC | vs degree baseline | dependence on # bound genes |
|---|---|---|---|
| Degree baseline (# genes the compound binds) | 0.701 (0.669-0.731) | - | 1.00 |
| DWPC, binding metapaths | 0.789 (0.755-0.820) | wins 83% of diseases | +0.43 |
| Permutation z-score, binding metapaths | 0.686 (0.646-0.723) | not different (p=0.63) | -0.48 |
| **Consensus** (mean of DWPC and z percentiles) | 0.758 (0.723-0.790) | wins 67% (p=0.001) | -0.06 |

* The graph score beats a "promiscuous drugs have more indications" baseline, but a large part of its accuracy
  is degree-related; the consensus keeps most of the accuracy with no degree dependence.
* Expression-signature evidence did not help (AUROC 0.64) and is used only as an annotation.
* Candidate lists are audited by **gene ablation**: some top candidates collapse when one gene is removed
  (e.g. rank 1 -> 162), others do not. See the evidence cards and `fig_ablation.png`.
* Binding edges carry **no direction of effect**: agonists and antagonists, toxicants and drugs look alike.

Figures: `outputs/figures/` (benchmark, damping sweep, hub audit, ablation, evidence cards, path graphs).
Interactive path explorer: `docs/interactive/index.html` (static page; click a candidate and a node to see its paths).
Literature triage of the shortlisted candidates (trial history, direction of effect): [`docs/literature_triage.md`](docs/literature_triage.md).

## Method in one paragraph

Hetionet edges are loaded into sparse matrices. For each of three compound-to-disease metapaths (compound binds gene
associated with disease; ... interacting with a disease gene; ... sharing a pathway with a disease gene) the
degree-weighted path count (damping w = 0.6) is computed for all 1,552 x 137 pairs. A z-score against 100
degree-preserving edge-swap permutations measures specificity. Percentile ranks of the two are averaged into a
consensus. Supporting paths for the shortlist are retrieved with **SPARQL (rdflib)** from an RDF export of the 13
mechanistic relations (362,912 triples) and their weights are checked against the matrix calculation.

## Reproduce

```bash
pip install -r requirements.txt
# put hetionet-v1.0.0.zip (https://zenodo.org/records/268568) in data/raw/
python3 src/run_all.py            # ~10 min on a laptop; --n-perm 20 for a quick run
```

| Step | Script | Output |
|---|---|---|
| 1-2 | `01_download.py`, `02_load_hetionet.py` | integer-indexed parquet (2.25M edges in ~1.6 s) |
| 3 | `03_schema_report.py`, `03b_rdflib_timing.py` | schema tables/figures, rdflib scale test |
| 4-5 | `04_build_rdf.py`, `05_sparql_queries.py` | RDF export, SPARQL queries |
| 6 | `06_metapath_scoring.py` (+ `dwpc.py`) | DWPC matrices and rankings |
| 7 | `07_evaluate.py`, `07b_pathway_leakage_check.py` | benchmark, baselines, damping, robustness |
| 8-9 | `08_permutation_null.py`, `09_consensus_candidates.py` | z-scores, consensus candidates |
| 10 | `10_explain_paths.py` (+ `data/triage.csv`) | paths, contribution shares, ablation |
| 11-12 | `11_figures.py`, `12_interactive_demo.py` | figures, static interactive explorer |

## Data and license

Hetionet's content is CC0, but it integrates many sources with their own licenses (some non-commercial); see the
Hetionet repository. The dataset is not redistributed here - the scripts read the Zenodo bundle you download.
`data/triage.csv` is a hand-made classification of compounds (endogenous, cytotoxic, direction, trial history); notes marked in `docs/literature_triage.md` were checked against searched sources, the rest come from general pharmacology.

## Limitations (short)

Hetionet is a 2017 snapshot built from 2015-16 sources; binding has no direction of effect; negatives in the
benchmark are unlabeled; Alzheimer's has 4 labeled "treats" drugs and Parkinson's none (28 "palliates"), so
disease-specific recovery is descriptive; w = 0.6 was chosen using the same benchmark. Details in the write-up.

## Citation

Himmelstein DS et al. *Systematic integration of biomedical knowledge prioritizes drugs for repurposing.*
eLife 2017;6:e26726.
