# Explainable knowledge-graph drug repurposing for Alzheimer's and Parkinson's disease

## 1. Motivation

Drug repurposing looks for new indications among compounds whose safety is already known. Biomedical knowledge
graphs make the search systematic: a compound and a disease are linked through genes, pathways and other diseases,
and the strength of those links can be scored. Two things make such scores hard to use. They are often reported as
a ranked list of names with no reason attached, and they can reflect properties of the graph (how well studied a
compound is) rather than biology. This project asks two questions of Hetionet (Himmelstein et al., 2017): does a
path-based score add anything beyond a simple popularity baseline, and can each candidate be shown as the
concrete evidence behind it?

In the Hetionet snapshot (2017), the drugs labeled as treating Alzheimer's disease are four symptomatic agents
(donepezil, galantamine, memantine, rivastigmine); Parkinson's disease has no "treats" edges at all and 28
"palliates" edges. Recovering these labels therefore says little about finding disease-modifying candidates, and
the evaluation is built around that limitation.

## 2. Data

Hetionet v1.0: 47,031 nodes of 11 types and 2,250,197 edges of 24 types (137 diseases, 1,552 compounds, 20,945
genes, 1,822 pathways). Files were loaded from the Zenodo release (edge list 11.8 MB compressed). The full graph
loads into pandas in under two seconds. An rdflib test showed the full graph is feasible in memory (~2 GB), but
unanchored SPARQL joins scale badly, so the RDF layer is restricted to the 13 mechanistic relations (320,145
edges, 362,912 triples), which drop the large annotation relations (gene-process, anatomy, side effects, ...).

## 3. Method

**Metapaths.** Three compound-to-disease paths: compound-binds-gene-associated-with-disease (CbGaD); binds a gene
that interacts with a disease-associated gene (CbGiGaD); binds a gene that shares a pathway with a
disease-associated gene (CbGpPWpGaD; paths that revisit the same gene are subtracted). An expression-signature
metapath (reversal minus concordance) and two similarity-to-known-treatment metapaths were also implemented.

**Degree-weighted path count (DWPC).** Each edge is weighted by (source degree x target degree)^-w within its
relation; a path's weight is the product along the path and a metapath score is the sum. w = 0.6, so that
generic hub genes and promiscuous compounds contribute less. Scores for all 1,552 x 137 pairs are computed with
sparse matrix products in under a second.

**Specificity (z-score).** 100 degree-preserving edge-swap permutations of the eight relations used by the
mechanistic metapaths give a null mean and standard deviation for every compound-disease pair; z = (observed -
null mean) / sqrt(null sd^2 + floor^2), with a variance floor for pairs that are almost always zero in the null.

**Aggregation and consensus.** Scores become coverage-aware percentile ranks (compounds without data for a
metapath are not ranked on it); DWPC and z aggregates are averaged into a consensus.

**Explanations.** For a shortlist, every supporting path is retrieved by SPARQL (rdflib) and re-weighted; the
path weights sum exactly to the matrix score for all shortlisted candidates. Contribution shares identify the
driving genes and pathways, and a counterfactual removes the top one or three disease genes and recomputes the
candidate's rank.

## 4. Evaluation

**Benchmark.** For each of 58 diseases with at least 3 "treats" edges (731 positives), all compounds are ranked
and AUROC computed; results are averaged over diseases with a bootstrap over diseases. Negatives are unlabeled
(some are true indications), so AUROC is conservative. A positive's own edge is never used to score it.

| Method | AUROC (95% CI) | Wins vs degree | Spearman with degree |
|---|---|---|---|
| Degree baseline | 0.701 (0.669-0.731) | - | 1.00 |
| DWPC, binding metapaths | 0.789 (0.755-0.820) | 83% | +0.43 |
| z-score, binding metapaths | 0.686 (0.646-0.723) | 47% (p = 0.63) | -0.48 |
| Consensus (mean) | 0.758 (0.723-0.790) | 67% (p = 0.001) | -0.06 |
| Consensus (min) | 0.742 (0.707-0.775) | 64% (p = 0.02) | -0.04 |

Single metapaths: CbGaD 0.720 (not better than degree, p = 0.34), CbGiGaD 0.691, CbGpPWpGaD 0.763, signature
0.636 (below degree). Damping barely changes accuracy (0.771 at w = 0 to 0.783 at w = 0.6-0.8 for the aggregate)
but lowers dependence on degree from 0.65 to 0.34.

**Interpretation.** Known indications correlate with compound degree (the degree baseline alone reaches 0.70),
so part of the DWPC advantage is popularity. Removing degree entirely (z-score) leaves a signal about as
predictive as degree but independent of it (negative correlation), and the consensus keeps most of the DWPC
accuracy with no degree dependence. The consensus does not remove hubs from the extreme top of the lists: for
Parkinson's the top-20 candidates bind a median of 25 genes versus 5 for a typical compound (Alzheimer's: 9).

**Robustness.** Dropping the 15 pathways named after a disease (e.g. "Alzheimers Disease") changes nothing
(AUROC 0.789 either way, 14/15 top candidates unchanged). Disease-named pathways are hidden from the path
figures because they are near-tautological as explanations.

**Alzheimer's and Parkinson's.** Treats-or-palliates drugs (9 and 28) have a median consensus rank percentile of
0.27 and 0.17 (random 0.5). This is descriptive only, and the known drugs are symptomatic agents acting on
targets already in the disease gene sets.

## 5. What the explanations show

Gene ablation separates candidates whose score rests on one gene from those that do not:

| Disease | Candidate | DWPC rank | without top gene | without top 3 genes |
|---|---|---|---|---|
| AD | Topiramate | 2 | 13 | 60 |
| AD | Marimastat | 3 | 261 | 293 |
| AD | Varenicline | 4 | 20 | 537 |
| AD | Minocycline | 5 | 17 | 274 |
| AD | Lithium | 7 | 8 | 321 |
| AD | Nicotine | 8 | 15 | 68 |
| AD | Clenbuterol | 10 | 65 | 110 |
| AD | Tyloxapol | 11 | 307 | 902 |
| PD | Dronedarone | 1 | 162 | 205 |
| PD | Sunitinib | 3 | 6 | 19 |
| PD | Topiramate | 4 | 18 | 411 |
| PD | Minocycline | 5 | 11 | 21 |
| PD | Vitamin E | 7 | 11 | 22 |
| PD | Metformin | 8 | 12 | 378 |
| PD | Sucralfate | 9 | 438 | 696 |
| PD | Isoprenaline | 11 | 329 | 344 |

Ranks are among 1,389 compounds with binding data (DWPC on the binding metapaths, w = 0.6). Bacitracin, the top
Alzheimer's consensus hit, is excluded from this shortlist by the triage step (probable wrong direction; see below).

The top Parkinson's candidate, dronedarone, gets 93% of its direct-binding score from a single gene (HCN3) and drops
to rank 162 without it. Sunitinib, in contrast, is supported by 132 bound genes and 6,870 pathway paths; it is
stable under ablation, but that is breadth from a promiscuous kinase inhibitor and not specific evidence. A
robust rank is therefore not the same as a specific one. Recurring driver genes group candidates into mechanistic
clusters (for Alzheimer's, nicotinic receptor genes; for Parkinson's, monoamine-system genes).

**Literature triage.** The shortlisted candidates were checked against published trials and reviews
(`docs/literature_triage.md`). Several of the best-supported graph hypotheses had already been tested clinically:
lithium (mixed), minocycline (negative in mild Alzheimer's; inconclusive futility result in early Parkinson's), nicotine
(positive pilot in mild cognitive impairment, larger trial registered) and varenicline (negative). Two top hits
appear to act in the wrong direction: bacitracin is an experimental inhibitor of an amyloid-degrading enzyme, and
tetrabenazine-type monoamine depleters are warned to worsen parkinsonism. The graph therefore mostly rediscovers
existing hypotheses and cannot by itself say which of them will work.

## 6. Limitations

* **No direction of effect.** "Binds" does not distinguish agonists from antagonists or drugs from toxicants.
  The lists contain examples where the mechanism likely runs opposite to benefit (for instance monoamine-depleting
  reserpine-class agents ranked for Parkinson's, and neuromuscular nicotinic blockers ranked beside nicotinic
  agonists for Alzheimer's). Direction would need curated pharmacology.
* **Non-therapeutic entries.** Endogenous cofactors (ATP, glutathione, FAD, choline), a diagnostic tracer, and
  cytotoxic or metal compounds appear among top hits; a hand-made triage table
  (`data/triage.csv`) marks them and they are excluded from the explained shortlist.
* **Snapshot and source bias.** Hetionet dates from 2015-16 sources; well-studied genes and drugs are
  over-represented; the benchmark labels come from overlapping literature as the features (some circularity).
* **Small disease-specific evidence.** 4 treats drugs for Alzheimer's, none for Parkinson's.
* **Tuning.** w = 0.6 was chosen on the benchmark used for evaluation (differences between w values were within
  noise); the z-scores use a variance floor and a finite null.
* **No external validation.** The benchmark uses Hetionet's own labels. A time-split test on trials that started after
  the graph's cut-off (about 2016) would be the natural next check; it was not done here. A source with action types
  and contraindications (for example PrimeKG) would address the direction-of-effect gap.
* **Similarity metapaths** (compound resembles a treating compound; treats a resembling disease) score well on the
  benchmark (the latter reaches 0.94) because the disease list contains near-synonymous diseases, and they are
  unavailable for the two target diseases. They are excluded from the headline results.

## 7. Relation to experimental validation

A graph score prioritizes hypotheses; it does not test them. Graph-derived support is literature-derived support
and mostly restates known pharmacology. Turning a candidate into a finding requires: confirming the direction of
effect at the implicated targets; mechanistic assays in disease-relevant cells; efficacy in animal models;
pharmacokinetics (brain exposure); and clinical evidence, ideally with retrospective analysis of real-world data.
The graph's strongest contribution is to make the evidence for each hypothesis inspectable, so that fragile or
tautological support can be filtered before expensive validation.

## 8. Reference

Himmelstein DS, et al. Systematic integration of biomedical knowledge prioritizes drugs for repurposing.
eLife 2017;6:e26726. Hetionet: https://github.com/hetio/hetionet
