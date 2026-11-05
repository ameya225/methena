# What the pipeline found in the schizophrenia datasets

This summary covers the public GEO data from Kiltschewskij, Reay & Cairns (2024, *Mol Psychiatry*,
doi:10.1038/s41380-024-02749-5), analysed with methena 1.0.0. Every number below comes from
the result tables under `results/<dataset>/analysis/`. Per-dataset details are in each
`analysis_report.html`.

## Datasets processed

| Dataset | Tissue / platform / input | Samples after QC (SCZ / control) | CpGs tested |
|---|---|---|---|
| GSE84727 (UCL) | whole blood, 450K, M/U intensities + detP | 819 (399 / 420); 28 PCA outliers removed | 437,075 |
| GSE89702 | cerebellum, 450K | 31 | 416,164 |
| GSE89703 | hippocampus, 450K | 26 (14 / 12) | 415,531 |
| GSE89705 | striatum (DBCBB), 450K | 29 | 416,043 |
| GSE89706 | striatum (LNDBB), 450K | 48 | 418,898 |
| GSE74193 subset | DLPFC, 450K **IDAT** (30 + 30 adults drawn) | 51 (27 / 24) | 434,125 |
| GSE166255 subset | whole-blood **RRBS** (neuropathic pain vs control; not SCZ; tests the sequencing input) | 24 (12 / 12) | 1,052,568 |

## 1. Blood (GSE84727): case-control differences

Model: M-value ~ group + DNAm age + sex + 6 blood-cell proportions (EpiDISH; neutrophils as the reference)
+ smoking score + slide + array row. 176 samples (150 of them cases) have no reported age, so we
followed the paper and used methylation-derived age (Zhang 2019 clock; MAE 2.6 years against reported age).

* **Inflation.** Before correction λ = 2.19 and 29,008 CpGs pass FDR < 0.05. The leading principal
  component is associated with diagnosis (P = 2×10⁻²²) and with lymphocyte proportions, so blood
  composition differs between groups beyond what the 7-cell reference captures. After bacon
  correction (estimated inflation 1.38), **4,294 DMPs** remain at FDR < 0.05, together with **505 DMRs**.
* **Direction and size.** 91 % of DMPs are hypermethylated in schizophrenia. Effects are small
  (median |Δβ| = 1.2 percentage points), as expected for a blood EWAS of a psychiatric disorder.
* **Where the changes are.** Hypermethylated DMPs are depleted from CpG islands (OR 0.05) and
  promoters (OR 0.47) and enriched in open sea (OR 3.7) and **enhancers (OR 2.5)**.
* **Transcription factors.** The sequences around hypermethylated DMPs are enriched for **AP-1
  (FOS/JUN/JUNB/JUND/FOSL1/2), BATF and SMAD2/3 motifs** compared with GC-matched background
  (monaLisa, JASPAR2020). AP-1/BATF enhancers are activity-dependent regulatory elements in
  myeloid and activated lymphoid cells. Methylation gain at these sites is expected to reduce
  TF binding and enhancer activity, which links the DMPs to immune and inflammatory gene regulation.
* **Regions.** The strongest DMRs lie in the MHC class III region on chr6, including
  ***TNF/LTA*** (two DMRs) and *HSD17B8/SLC39A7/MIR219A1*, plus *WBP1L* and *RPL13A*. These
  results are consistent with the established MHC and immune involvement in schizophrenia.
* **Genes and predicted transcriptional effects.** 900 genes carry promoter hypermethylation
  (predicted lower expression, e.g. *GPR37L1, MIR134, MIR382, FGF1, RUNX1T1*). 127 genes carry promoter
  hypomethylation (predicted higher expression, e.g. *CORO1C, PDE7B, VGLL4, KLF3*).
  These are predictions; `analysis.expression_results` can test them against paired RNA-seq.
* **Pathways** (missMethyl, corrected for CpGs-per-gene bias). GO cell-periphery, plasma-membrane,
  receptor-signalling and morphogenesis terms. Reactome *Signaling by GPCR* and *G-protein-gated
  potassium channels* (FDR 0.028). KEGG *Rap1 signalling*, *aldosterone synthesis* and *dopaminergic
  synapse* (FDR 0.010-0.017).
* **Variability (the focus of the paper).** Only **3 VMPs** reach FDR < 0.05 in this single
  cohort (*ZFPM1, SPTBN5, MITF*; all more variable in cases), and none of them is a DMP. The paper's 213
  VMPs came from a three-cohort meta-analysis of about 2,000 people, so one cohort is expected to
  detect few. The pipeline does reproduce the qualitative conclusion that **variance and mean effects
  mark different CpGs**.

## 2. Epigenetic clocks (blood, 643 samples with reported age)

* **Benchmark.** Zhang 2019 (r = 0.98, MAE 2.6 y) and skin-and-blood Horvath 2018 (r = 0.97,
  MAE 2.6 y) are the best calibrated. Horvath 2013 MAE is 4.6 y, Hannum 7.4 y (systematic +7 y bias),
  and PC-GrimAge correlates well (r = 0.93) but has a +16.6 y offset (as expected; it is calibrated as a mortality-risk age).
  The cortical clock (Shireby) is the best clock in all three brain regions (hippocampus MAE 3.4 y,
  striatum 4.0-4.5 y). This tissue specificity of clock accuracy confirms that the benchmark works.
* **Case-control.** After adjusting for sex, cells, smoking and batch, schizophrenia is associated with
  **faster pace of ageing (DunedinPACE +0.039 years/year, +0.28 SD, FDR 3×10⁻⁴)**, higher
  **PC-PhenoAge (+0.82 y, FDR 0.026)** and skin-and-blood age (+0.65 y, FDR 0.033). It is also associated with lower
  mitotic-clock scores (epiTOC/epiTOC2, -0.21 SD, FDR 0.012). First-generation chronological
  clocks (Horvath 2013, Hannum) do not differ. This is the expected pattern when a disease
  affects *biological* (morbidity-linked) ageing but not the chronological-age component.

## 3. Brain and cross-tissue concordance

* With 26-48 donors per region, no CpG reaches FDR < 0.05 in any brain region (λ 0.94-1.02,
  well calibrated). A handful of VMPs are found (0-6 per region).
* **Blood DMPs replicate weakly in brain.** For the 4,294 blood DMPs, the sign of the case-control
  effect agrees in **54 % in hippocampus (binomial P = 2×10⁻⁶)** and **54 % in LNDBB striatum
  (P = 3×10⁻⁸)**, but not in cerebellum (49 %) or DBCBB striatum (51 %). Only 114 (hippocampus) and 103
  (LNDBB striatum) blood DMPs replicate at P < 0.05 in the same direction, against about 105 expected
  by chance. The shared signal is therefore a slight overall shift in direction in two of four regions,
  not a set of individually replicating CpGs, which is typical for blood-brain comparisons. This mirrors the paper's report of concordance between blood VMPs and some brain regions.
* **DLPFC (IDAT path).** Without correction 10,241 CpGs pass FDR < 0.05, but λ = 2.22 and bacon
  removes all of them (minimum corrected P = 7.5×10⁻⁷). The PC diagnostics show that PC1 is the
  neuronal proportion and PC4 tracks reported race. In a 51-sample subset this pattern most likely reflects residual
  confounding, not disease. The report flags it, and `n_pcs` or batch covariates (e.g. plate) are
  the recommended next step. Blood DMPs do **not** agree in direction in DLPFC (45.5 % concordant, i.e. slightly discordant; binomial P = 0.005).
  Gene-set signals in DLPFC are nominal (GABA-ergic synapse, axonogenesis; FDR > 0.05).

## 4. Bisulfite-sequencing example (GSE166255 RRBS)

The sequencing path (Bismark coverage → coverage filters → count-based M-values with coverage
weights → limma/DMRcate/varFit) ran end to end on 1.05 M CpGs. Four of the 12 cases are profile
outliers. Removing them would have taken 33 % of one group, so the safeguard **kept and flagged**
them. No DMP passes FDR < 0.05, and the 484 VMPs (λ = 2.35) are probably driven by those outliers.
The report says so explicitly.

## Caveats

* Associations are cross-sectional. Antipsychotic treatment (clozapine affects neutrophils), smoking
  and cell composition are strong confounders in blood. The adjustment reduces but does not remove them,
  as the inflation shows. All results should be treated as associations, not causes.
* Promoter-methylation predictions of transcription and TF-motif enrichment are hypotheses that
  need expression or chromatin data to confirm.
* The GSE74193 and GSE166255 runs use subsets chosen to test input types, not to make biological claims.
