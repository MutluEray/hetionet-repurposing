# Literature triage of the shortlisted candidates

Searches run on 29 Sep 2026 (web search; abstracts, registries and review articles). This is a **non-systematic
check of the specific claims the write-up depends on**, not a literature review. "No literature found" means
nothing turned up in these searches, not that none exists. Everything below should be re-checked against the
primary papers before it is quoted elsewhere.

Categories match `data/triage.csv`. "Graph evidence" is what our pipeline reports (`explain_summary_*.csv`).

## Alzheimer's disease

| Candidate | Graph evidence | Literature status | Assessment |
|---|---|---|---|
| Bacitracin | binds IDE and A2M; 2 bound genes | Used experimentally as an inhibitor of insulin-degrading enzyme (IDE); IDE degrades amyloid-beta ([GeneCards IDE](https://www.genecards.org/card/IDE); Qiu et al., J Biol Chem 1998;273:32730) | **wrong direction likely**: inhibiting an amyloid-degrading enzyme runs opposite to benefit. Excluded from the explained shortlist |
| Marimastat | 100% of direct-binding score from ADAM10 | ADAM10 is the major alpha-secretase and increasing its activity is the proposed strategy ([Frontiers Mol Neurosci 2017](https://www.frontiersin.org/journals/molecular-neuroscience/articles/10.3389/fnmol.2017.00056/pdf)); hydroxamate metalloproteinase inhibitors such as batimastat block alpha-secretase (Hooper & Turner, Curr Med Chem 2002;9:1107) | **direction unclear, probably opposite**; marimastat's own ADAM10 activity was not verified in these searches |
| Lithium | GSK3A, generic pathways; 5 bound genes | 10-week RCT in mild AD (n=71) found no effect on GSK-3 activity or CSF biomarkers ([Hampel et al., J Clin Psychiatry 2009;70:922](https://www.psychiatrist.com/jcp/lithium-trial-alzheimers-disease-randomized-single/)); MCI trials by Forlenza et al. (Br J Psychiatry 2011, 2019) exist; overall mixed | **investigated, mixed**; mechanistically coherent (GSK-3) |
| Minocycline | CASP3, CYCS, generic apoptosis genes | MADE phase II RCT, 554 mild-AD patients, 24 months: did not slow cognitive or functional decline (Howard et al., [JAMA Neurol 2020;77:164](https://doi.org/10.1001/jamaneurol.2019.3762)) | **investigated, negative** |
| Nicotine | CHRNB2, CHRNA4 (nicotinic receptor genes) | 6-month MCI pilot RCT (n=74): improved attention/cognitive tests, no significant clinician-rated global effect ([Newhouse et al., Neurology 2012;78:91](https://scholars.duke.edu/publication/783728)); larger MIND trial registered ([NCT02720445](https://clinicaltrials.gov/study/NCT02720445)); its results were not found | **investigated, pilot signal** |
| Varenicline | CHRNB2 | Phase II crossover in mild-to-moderate AD: no improvement in cognition, behaviour or global change ([Kim et al., Dement Geriatr Cogn Disord 2014;37:232](https://kct.medric.or.kr/Controls/PVIEW.aspx?i=24247022)); an [Alzforum summary](https://alzforum.org/node/2787) reports worse neuropsychiatric state and GI effects | **investigated, negative** |
| Topiramate | kainate/AMPA-receptor genes (GRIK1/2, GRIA1) | No AD trial found. Cognitive adverse effects, especially verbal fluency and word-finding, are well documented ([Sommer et al., Ther Adv Neurol Disord 2013;6:211](https://pmc.ncbi.nlm.nih.gov/articles/PMC3707352)) | **tolerability concern** in a cognitive disease |
| Tyloxapol | binds one gene (LPL) | Nothing relevant found | **no evidence found**; single-gene, fragile |
| Cytisine | nicotinic receptor genes | No AD trial found | no evidence found |
| Vitamin E | GSTA4, oxidative-stress genes | Slowed progression in moderately severe AD (Sano et al., NEJM 1997;336:1216); no benefit for MCI conversion (Petersen et al., NEJM 2005;352:2379); TEAM-AD (n=613): slower functional decline, no cognitive benefit ([Dysken et al., JAMA 2014;311:33](https://sites.bu.edu/geriatricsfellowship/files/2014/06/02.14.14_Lee_-Vit-E-for-AD.pdf)) | **investigated, mixed** |

## Parkinson's disease

| Candidate | Graph evidence | Literature status | Assessment |
|---|---|---|---|
| Minocycline | CYCS, CASP3, IL1B | NET-PD futility trial (n=200): could not be rejected as futile, i.e. eligible for further study, not evidence of efficacy; tolerability was lower and 23% discontinued by 18 months ([Neurology 2006;66:664](https://scholars.duke.edu/publication/1175138); [Clin Neuropharmacol 2008;31:141](https://pubmed.ncbi.nlm.nih.gov/18520981/)) | **investigated, inconclusive** |
| Sunitinib | 132 bound genes; kinome-wide overlap | Nothing found for sunitinib. The class hypothesis (c-Abl inhibition) was tested with nilotinib: protection in an MPTP mouse model ([Karuppagounder et al., Sci Rep 2014;4:4874](https://pmc.ncbi.nlm.nih.gov/articles/PMC4007078)); a small uncontrolled safety trial whose efficacy signals were questioned ([J Parkinsons Dis commentary](https://journals.sagepub.com/doi/abs/10.3233/JPD-160904)); a meta-analysis of 3 RCTs (163 patients) found no significant difference from placebo on tolerability, adverse events or HVA ([Front Aging Neurosci 2022](https://www.frontiersin.org/journals/aging-neuroscience/articles/10.3389/fnagi.2022.996217/epub)) | class hypothesis is live but weakly supported; the graph's support for sunitinib is **broad, not specific** |
| Metformin | ND3, oxidative phosphorylation | Observational only and inconsistent: one meta-analysis found no effect on neurodegenerative disease overall and higher PD risk with monotherapy (OR 1.66; [Ping et al. 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7390234)); another found lower dementia risk (RR 0.79) but no PD/AD reduction ([Zhang et al. 2022](https://doi.org/10.1177/20406223221109454)); [Alzforum review](https://alzforum.org/node/1324341); MAP RCT in aMCI registered ([NCT04098666](https://clinicaltrials.gov/study/NCT04098666)) | **investigated, confounded and mixed** |
| Vitamin E | GSTA4, dopamine metabolism | Its DATATOP arm served as the placebo-like calibration group when NET-PD set its futility threshold; no efficacy claim | investigated (AD literature above) |
| Dronedarone | 93% of direct score from HCN3 | Blocks human HCN4 channels in vitro, IC50 about 1 uM ([Naunyn Schmiedebergs Arch Pharmacol 2011;383:347](https://www.medchemexpress.com/publications/21279331.html)); no neuro literature found | **no evidence found**; single-gene, fragile |
| Isoprenaline, Sucralfate | ADRBK1; FGB (72% of direct score) | Nothing relevant found | no evidence found; single-gene, fragile |
| Tetrabenazine (earlier top-25) | VMAT/monoamine genes | VMAT2 inhibitor that depletes monoamines; product monographs state it can induce parkinsonism and exacerbates parkinsonian symptoms ([monograph](https://pdf.hres.ca/dpd_pm/00057879.PDF)); central effects resemble reserpine ([monograph](https://pdf.hres.ca/dpd_pm/00057874.PDF)) | **wrong direction** |
| Deserpidine, Rescinnamine, Metyrosine | monoamine-system genes | Reserpine-class and tyrosine-hydroxylase-inhibiting agents lower monoamines; class membership is from general pharmacology and was not verified in these searches | **wrong direction likely** |

## What this changes

* Bacitracin (rank 1 in Alzheimer's) is excluded from the explained shortlist as a probable wrong-direction hit.
* The two best-supported "graph recovers a real hypothesis" cases (lithium, nicotinic agonists) were already
  tested clinically with mixed or negative results; minocycline and varenicline were negative in Alzheimer's.
  The graph is therefore mostly rediscovering known hypotheses, which is what a graph built from the literature
  should do, and it does not identify winners.
* Repurposing hypotheses from the graph need the direction of effect (agonist vs antagonist, inhibitor vs
  activator) before they mean anything; this is the largest gap.
