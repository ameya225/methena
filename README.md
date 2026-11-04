# methylflow: a reproducible DNA-methylation pipeline

[![R-CMD-check](https://github.com/ameya225/methylflow/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ameya225/methylflow/actions/workflows/R-CMD-check.yaml)

methylflow takes raw DNA-methylation data and turns it into a report that says
**where and how cases differ from controls**, which genes and pathways are involved,
how this could change gene transcription, and whether cases look biologically older.
You run it with one command, and the same input with the same settings gives the
same result tables every time.

```
 raw data ──► 1. PREPROCESS ─────────────────────► 2. ANALYZE ─────────────────────────► 3. EXPLORE
 (IDAT, intensity    QC, filtering, normalisation,     differential methylation (CpGs,     HTML report +
  matrices, betas,   cell-type / sex / smoking         regions, variability), 17 epigenetic interactive
  or Bismark files)  estimates  ► methylflow_dataset   clocks, genes, pathways, TF motifs   Shiny app
```

It supports:

| Your data | `input.type` | Typical source |
|---|---|---|
| Illumina IDAT files (450K, EPIC, EPICv2) | `idat` | the scanner, or GEO `*_RAW.tar` |
| Methylated + unmethylated signal tables (+ detection p-values) | `intensity` | GEO supplementary files, GenomeStudio |
| A table of beta values (or M-values) | `beta` | collaborators, GEO "processed" files |
| Bisulfite sequencing (WGBS, RRBS, EM-seq) coverage files | `bismark` | Bismark `*.bismark.cov.gz` |

---

## 1. What you need

* A Mac or Linux computer (Windows users: use WSL2/Ubuntu). 16 GB RAM is enough for up to
  ~300 450K samples; ~850 samples need 32 GB.
* **R 4.4 or newer**, from <https://cran.r-project.org> (on macOS use the official
  installer `.pkg`).
* About 10 GB of free disk for R packages and reference data, plus space for your data.
* Internet access during setup. Analyses then run offline, with three optional
  exceptions that are downloaded on first use: KEGG and MSigDB gene sets, and the JASPAR2024 motif
  database. If JASPAR2024 cannot be reached, the bundled JASPAR2020 motifs are used instead.

## 2. Install (once)

**Option A: download the whole pipeline (recommended for most users).** Open the **Terminal** app and run:

```bash
git clone https://github.com/ameya225/methylflow.git      # or: Code -> Download ZIP on GitHub, then unzip
cd methylflow
./run_methylflow.sh setup        # installs ~350 R packages into ./Rlib (30-60 min, only once)
./run_methylflow.sh demo         # 2-minute end-to-end check on simulated data
```

The demo should finish with `DEMO OK: recovered 9x% of the 200 planted differences`
and print the path of an HTML report you can open in a browser.
The packages go into the `Rlib` folder inside this directory, so the R installation
you use for other work is not modified.

**Option B: install only the R package (for R users).**

```r
install.packages(c("remotes", "BiocManager"))
source("https://raw.githubusercontent.com/ameya225/methylflow/main/inst/scripts/install_dependencies.R")
remotes::install_github("ameya225/methylflow")
library(methylflow)
```

The command-line script is then at `system.file("scripts", "methylflow.R", package = "methylflow")`.

## 3. Prepare your data

1. Put your raw files in one folder (e.g. `my_data/`).
2. Make a **sample sheet**: a CSV or tab-separated file with one row per sample. Minimum columns:

   | sample_id | group | age | sex | (data id column) |
   |---|---|---|---|---|
   | P001 | case | 54 | F | 200123450001_R01C01 |
   | P002 | control | 61 | M | 200123450001_R02C01 |

   * `group` must have exactly two values (for example case/control).
   * `age` (years) is needed for the epigenetic-clock comparison; `sex` (M/F) is used for
     adjustment and for the automatic sex check.
   * The data id column tells methylflow which column of the data file, IDAT file
     name or coverage file belongs to the sample. It is detected automatically when
     possible. IDATs named like `GSM123_2001234_R01C01_Grn.idat.gz` are matched
     to the sample id automatically.
   * Data from GEO: instead of a sample sheet you can point `input.geo_series_matrix` to the
     `GSExxxx_series_matrix.txt.gz` file; every "characteristics" field becomes a column.
     To see these columns as a table, run
     `Rscript inst/scripts/methylflow.R geo-sheet --series GSE..._series_matrix.txt.gz --out samples.tsv`.

## 4. Write a configuration file

Copy the template and edit the lines marked `EDIT`:

```bash
cp inst/config/config_template.yaml my_study.yaml
```

Every option is explained in the template. Paths are relative to the YAML file. The
`configs/` folder contains working examples for each input type (see section 8).

Check the configuration before a long run:

```bash
./run_methylflow.sh check my_study.yaml
```

This reads the sample sheet, applies your group/age/sex settings and prints how many
samples are in each group, without processing any data.

## 5. Run

```bash
./run_methylflow.sh all my_study.yaml          # = preprocess + analyze
# or step by step
./run_methylflow.sh preprocess my_study.yaml
./run_methylflow.sh analyze    my_study.yaml   # re-run with new settings without re-doing QC
```

Measured run times on an 8-core laptop with 32 GB RAM:

| Dataset | preprocess | analyze |
|---|---|---|
| 26-48 brain samples (450K intensities) | ~1 min | ~1.5-2 min |
| 60 IDATs (450K) | ~1 min | ~4 min |
| 24 RRBS samples (1 M CpGs) | ~2 min | ~3 min |
| 847 blood samples (450K intensities) | ~18 min | ~40-75 min (the 17 clocks take most of this; set `threads:` or shorten `clock_set`) |

A log file is written next to the results.

## 6. Look at the results

```
<output_dir>/
├── preprocess/
│   ├── qc_report.html                 ◄ open first: sample QC, excluded samples, PCA, cell types
│   ├── sample_qc.tsv                  per-sample metrics and exclusion reasons
│   ├── probe_filter_summary.tsv       how many CpGs each filter removed
│   ├── sample_annotation_final.tsv    phenotype + estimated cell types, smoking score, PCs
│   ├── methylflow_dataset.rds         standardised data object (input to the analysis step)
│   └── run_manifest.json              inputs/outputs with SHA-256 checksums, settings, package versions
└── analysis/
    ├── analysis_report.html           ◄ the main report, starting with "Key findings"
    ├── dmp_results.tsv.gz             every CpG: effect size, P, FDR, gene, island context
    ├── dmr_results.tsv                differentially methylated regions
    ├── vmp_results.tsv.gz             differentially *variable* CpGs
    ├── gene_level_regulation.tsv      per gene: promoter vs body change and predicted effect on transcription
    ├── enrichment_dmp_GO/KEGG/Reactome.tsv   pathways (corrected for CpGs-per-gene bias)
    ├── regulatory_context_enrichment.tsv     island / promoter / enhancer enrichment of DMPs
    ├── tf_motif_enrichment.tsv        transcription-factor motifs near DMPs (if motif: true)
    ├── clock_estimates.tsv            17 epigenetic clocks per sample
    ├── clock_benchmark.tsv            accuracy of each clock vs chronological age
    ├── clock_group_tests.tsv          age acceleration, cases vs controls
    ├── cross_dataset_concordance.tsv  replication in other datasets/tissues (if compare_with set)
    ├── pc_covariate_association.tsv   what drives the main axes of variation (confounding check)
    ├── results_bundle.rds             data for the Shiny app
    └── figures/                       all plots as PNG
```

### Interactive browser (Shiny)

```bash
./run_methylflow.sh shiny <output_dir>/analysis
```

A browser tab opens with tabs for a volcano plot (click a point to see that CpG), a
Manhattan plot, filterable DMP/DMR/VMP/gene tables, a per-CpG viewer with values by group,
pathways, transcription-factor motifs, clocks and concordance. Use the sliders to change
the FDR and effect-size thresholds; **Download** saves the filtered table. Stop the app
with `Ctrl+C` in the Terminal.

## 7. How to read the main results

| Term | Meaning |
|---|---|
| **Beta value** | Fraction of DNA copies methylated at a CpG (0 = never, 1 = always). |
| **M-value** | log2(methylated/unmethylated). Used for statistics because it behaves better than beta. |
| **DMP** | Differentially methylated position: a CpG whose *average* methylation differs between groups. |
| **Δβ (delta beta)** | Difference in mean methylation, case minus control, adjusted for covariates. 0.02 = 2 percentage points. Blood EWAS of psychiatric traits typically find 1-3 %. |
| **DMR** | Differentially methylated region: several neighbouring CpGs changing together. These results are usually more robust than single CpGs. |
| **VMP** | Variably methylated position: methylation is more (or less) *variable* in cases. This can indicate heterogeneity among patients. |
| **FDR** | False discovery rate. At FDR < 0.05, about 5 % of reported hits are expected to be false. |
| **λ (lambda) / bacon** | Inflation of test statistics. λ ≈ 1 is ideal. λ > 1.2 indicates broad confounding (cell types, smoking, medication, batch) or a very polygenic signal. methylflow uses *bacon* by default to correct P-values for this, and reports both corrected and uncorrected counts. |
| **Promoter hypermethylation** | Usually silences the gene → "likely reduced transcription". Gene-body methylation is often *positively* linked to expression, so it is reported separately. |
| **Epigenetic age acceleration** | Clock age minus the age expected for the person's chronological age. Positive values mean the methylome looks older than expected. DunedinPACE is the *pace* of ageing (1.0 = one biological year per calendar year). |

The **Key findings** box at the top of `analysis_report.html` is generated from the results:
it states how many CpGs, regions and genes differ, which genomic contexts and pathways are
enriched, which genes are predicted to be repressed or de-repressed, which clocks differ
between groups and whether hits replicate in the comparison datasets.

## 8. The example study (schizophrenia, Kiltschewskij et al. 2024, *Mol Psychiatry*)

`configs/` reproduces analyses of the public GEO datasets used in that paper:

| Config | Dataset | Input path tested |
|---|---|---|
| `GSE84727_blood.yaml` | UCL blood cohort, 450K, 847 samples (414 SCZ / 433 controls) | intensity + detection p |
| `GSE89702_cerebellum.yaml`, `GSE89703_hippocampus.yaml`, `GSE89705_striatum.yaml`, `GSE89706_striatum.yaml` | post-mortem brain regions | intensity; GSM↔array linking by profile |
| `GSE74193_DLPFC_idat.yaml` | 60-sample adult DLPFC subset of GSE74193 | IDAT |
| `GSE166255_RRBS_blood.yaml` | whole-blood RRBS, 24-sample subset (not from the paper; tests the sequencing input) | Bismark coverage |

To re-create everything from scratch:

```bash
cd example_data
./download_geo.sh raw GSE84727:matrix GSE84727:GSE84727_detectionP.csv.gz \
   GSE84727:GSE84727_methylatedIntensities.csv.gz GSE84727:GSE84727_unmethylatedIntensities.csv.gz
# brain sets: see example_data/README_data.md for the full list of files
cd ..
./run_methylflow.sh all configs/GSE84727_blood.yaml
```

The brain intensity files are labelled with array IDs, while GEO metadata uses GSM accessions.
`inst/scripts/link_geo_samples.R` links the two by matching each sample's
methylation profile to the series-matrix values. In our runs all 142 brain samples
matched unambiguously (r = 1.00 against a runner-up of about 0.2-0.4).

## 9. Making sure it works (tests)

```bash
./run_methylflow.sh test
```

The automated tests check that:

* bad configs (missing files, wrong types, impossible thresholds) give clear error messages;
* simulated data with 150-200 planted differences are recovered (> 85 % sensitivity,
  < 10 % false discoveries) for array input and for sequencing input;
* QC removes a planted failing sample and a planted sex-mismatched sample;
* two runs give byte-identical result tables;
* beta and M-value inputs give the same dataset;
* missing covariates, missing values, a single group and mismatched sample sheets are handled;
* the memory-efficient dasen normalisation equals wateRmelon's implementation;
* the Shiny app builds and its filters work.

## 10. Troubleshooting

| Message | What to do |
|---|---|
| `Package 'X' is required ... run install_dependencies.R` | Run `./run_methylflow.sh setup` again; it only installs what is missing. |
| `Could not match any data columns to sample sheet rows` | Set `input.data_id_col` to the sample-sheet column that holds the names used in the data file (the message lists examples). |
| `group_col '...' not in sample sheet` | Check spelling; `./run_methylflow.sh check` prints all column names. |
| `Case-control analysis needs exactly 2 groups` | Set `phenotype.control_level`/`case_level`, or use `recode` to merge labels, or `keep_samples_query` to subset. |
| `vector memory limit ... reached` | Use a machine with more RAM, or analyse a subset (`keep_samples_query`). |
| Many samples excluded | Open `preprocess/qc_report.html`; thresholds are in the `qc:` section. |
| λ well above 1.2 | Look at the PC-covariate plot; add batch columns (`batch_cols`) or `n_pcs: 5`. |
| KEGG/MSigDB warnings | These gene sets need internet access. GO and local GMT files (`gene_set_files`) work offline. |
| `JASPAR2024 database not reachable` | Harmless: the motif step uses the bundled JASPAR2020 motifs instead. |
| `... outliers - too many to remove automatically` | The samples are kept and flagged. Look at them in the QC report; set `qc.force_drop_outliers: true` if they really are failed samples. |
| `strands not collapsed` (sequencing) | Install `BSgenome.Hsapiens.UCSC.hg38` (or hg19) to merge the two strands of each CpG. Results stay valid without it. |
| Any other error | Re-run with `METHYLFLOW_DEBUG=1 ./run_methylflow.sh ...` to print the call stack. |

## 11. What is in this repository

```
methylflow/
├── README.md               this file
├── run_methylflow.sh       one-command wrapper (setup / demo / check / all / shiny / test)
├── DESCRIPTION, NAMESPACE  R package metadata (the repository root is an installable R package)
├── R/                      package code (io, qc, normalize, covariates, dm, clocks, enrichment, analyze, shiny ...)
├── inst/scripts/           methylflow.R (command line), install_dependencies.R, link_geo_samples.R
├── inst/shiny/             Shiny app
├── inst/config/            config_template.yaml (every option documented)
├── inst/extdata/           bundled cross-reactive probe lists
├── tests/testthat/         automated tests (simulated data; no downloads needed)
├── man/                    function help pages
├── configs/                ready-to-run configs for the example GEO datasets
├── example_data/           download script, sample sheets and URL lists (raw data not included)
├── docs/                   INSIGHTS.md (example-study findings) and validation/ (test and check logs)
└── .github/workflows/      continuous integration: R CMD check (Bioconductor container) on every push
```

The package can also be used from R directly:

```r
library(methylflow)
cfg <- mf_read_config("my_study.yaml")
se  <- mf_preprocess(cfg)          # SummarizedExperiment: assay "beta", colData = phenotype + QC + covariates
res <- mf_analyze(cfg)             # list with $dmp, $dmr, $vmp, $clocks, $enrich, $reg, $motif, $genes
mf_shiny("results/my_study/analysis")
```

## 12. Methods and references

* QC/normalisation: minfi (Aryee 2014; Noob, Funnorm), wateRmelon dasen/BMIQ (Pidsley 2013; Teschendorff 2013);
  cross-reactive probes (Chen 2013; McCartney 2016; Pidsley 2016); dbSNP common SNPs from the Illumina annotation.
* Cell composition: EpiDISH (Teschendorff 2017; blood), HiBED (Zhang 2023; brain).
* Statistics: limma (Ritchie 2015), bacon (van Iterson 2017), DMRcate (Peters 2015),
  missMethyl varFit/gometh/gsameth (Phipson 2014, 2016).
* Clocks (dnaMethyAge engine): Horvath 2013 and 2018, Hannum 2013, PhenoAge (Levine 2018),
  Zhang 2019, cortical clock (Shireby 2020), PedBE (McEwen 2019), PC clocks incl. PC-GrimAge
  (Higgins-Chen 2022), DunedinPACE (Belsky 2022), DNAmTL (Lu 2019), epiTOC/epiTOC2 (Yang 2016;
  Teschendorff 2020), mortality score (Zhang 2017).
* Motifs: monaLisa (Machlab 2022) with JASPAR CORE vertebrate PWMs (JASPAR2024 if reachable, otherwise the bundled JASPAR2020).
* Gene sets: GO (GO.db), KEGG (online), Reactome GMT, optional MSigDB (msigdbr).
