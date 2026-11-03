# methylflow validation log (2026-11-02)

## R CMD check (with tests)
* checking extension type ... Package
* checking examples ... NONE
* checking tests ...
Status: OK

## Real-data reproducibility: GSE89703 run twice from raw files (SHA-256, first 16 hex)
(Done with the build immediately before the final change, which applies bacon only when the estimated inflation is > 1. That change alters the values but not the determinism of the code.)
- preprocess/sample_qc.tsv  a8f3fb1adb296f28  a8f3fb1adb296f28  IDENTICAL
- preprocess/probe_filter_summary.tsv  6cae4921a22398d1  6cae4921a22398d1  IDENTICAL
- preprocess/sample_annotation_final.tsv  8977a2cb87781320  8977a2cb87781320  IDENTICAL
- analysis/dmp_results.tsv.gz  3c490d424a0e309b  3c490d424a0e309b  IDENTICAL
- analysis/vmp_results.tsv.gz  a1aa3ec0760a762c  a1aa3ec0760a762c  IDENTICAL
- analysis/gene_level_regulation.tsv  1bf8be1b18d79bb5  1bf8be1b18d79bb5  IDENTICAL
- analysis/clock_estimates.tsv  93cda24e5a8c3f1a  93cda24e5a8c3f1a  IDENTICAL
- analysis/clock_benchmark.tsv  330ff6167f39c6fa  330ff6167f39c6fa  IDENTICAL
- analysis/regulatory_context_enrichment.tsv  12e3c258f553a812  12e3c258f553a812  IDENTICAL

## Blood GSE84727
Two independent analysis runs with the previous build gave identical counts (4,377 DMPs, 514 DMRs, 3 VMPs). The final build computes the t-to-z conversion exactly on the log scale; before, |t| > ~8 was capped. With the final build the counts are 4,294 DMPs, 505 DMRs and 3 VMPs.

## End-to-end runs (exit status 0 for all)
- GSE84727_blood: DMP 4294 | DMR 505 | VMP 3 | lambda 2.193
- GSE89702_cerebellum: DMP 0 | DMR 0 | VMP 6 | lambda 0.949
- GSE89703_hippocampus: DMP 0 | DMR 0 | VMP 3 | lambda 1.02
- GSE89705_striatum: DMP 0 | DMR 0 | VMP 0 | lambda 0.941
- GSE89706_striatum: DMP 0 | DMR 0 | VMP 3 | lambda 0.966
- GSE74193_DLPFC_idat: DMP 0 | DMR 0 | VMP 2865 | lambda 2.222
- GSE166255_rrbs: DMP 0 | DMR 0 | VMP 484 | lambda 1.258

## testthat output
```
array-pipeline: Setting options('download.file.method.GEOquery'='auto')
Setting options('GEOquery.inmemory.gpl'=FALSE)
.........................No methods found in package 'RSQLite' for request: 'dbListFields' when loading 'lumi'
...[2026-11-02 16:26:31] WARN  Covariate 'not_a_column' not available - skipped
.[2026-11-02 16:26:31] WARN  Covariate 'not_a_column' not available - skipped
.[2026-11-02 16:26:31] ERROR Case-control analysis needs exactly 2 groups; found: case
.[2026-11-02 16:26:31] ERROR No 'group' column: set phenotype$group_col in the config
.............[2026-11-02 16:27:03] INFO  Model: ~group + age_or_dnam + sex | n = 38 (control=18, case=20)
..
bisulfite: ....[2026-11-02 16:27:38] ERROR No coverage files matching '\.cov(\.gz)?$' in /Users/akulka77/.claude-science/orgs/39a431af-cd35-4f9c-a426-3289cd1d667b/workspaces/9d9d4f8b-d129-49c7-8e17-304b4d9ab64e/.rtmp/RtmpPx9TF3/file115b57c0a71b4/cov
.
config: [2026-11-02 16:27:38] ERROR Config file not found: /Users/akulka77/.claude-science/orgs/39a431af-cd35-4f9c-a426-3289cd1d667b/workspaces/9d9d4f8b-d129-49c7-8e17-304b4d9ab64e/.rtmp/RtmpPx9TF3/file115b56a3ba809/nope.yaml
.[2026-11-02 16:27:38] ERROR Config is not valid YAML: (/Users/akulka77/.claude-science/orgs/39a431af-cd35-4f9c-a426-3289cd1d667b/workspaces/9d9d4f8b-d129-49c7-8e17-304b4d9ab64e/.rtmp/RtmpPx9TF3/file115b56a3ba809/bad.yaml) Parser error: while parsing a flow sequence at line 1, column 8 did not find expected ',' or ']' at line 2, column 1
.[2026-11-02 16:27:38] ERROR input$type must be one of: idat, intensity, beta, bismark
.[2026-11-02 16:27:38] ERROR File/folder for input$sample_sheet does not exist: /Users/akulka77/.claude-science/orgs/39a431af-cd35-4f9c-a426-3289cd1d667b/workspaces/9d9d4f8b-d129-49c7-8e17-304b4d9ab64e/.rtmp/RtmpPx9TF3/file115b56a3ba809/missing.csv
.[2026-11-02 16:27:38] ERROR analysis$fdr must be in (0, 1)
.[2026-11-02 16:27:38] ERROR normalization$method must be one of: auto, noob, noob_bmiq, funnorm, quantile, swan, dasen, bmiq, none
...........
dasen: .

══ DONE ════════════════════════════════════════════════════════════════════════
🌈 Your tests are over the rainbow 🌈
```
