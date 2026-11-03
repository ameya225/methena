# Data used for the worked examples

All files are public (NCBI GEO). `download_geo.sh` downloads with resume/retry and logs SHA-256
checksums to `raw/MANIFEST.tsv`.

```bash
cd example_data
# UCL blood cohort (GSE84727)
./download_geo.sh raw GSE84727:matrix GSE84727:GSE84727_detectionP.csv.gz \
  GSE84727:GSE84727_methylatedIntensities.csv.gz GSE84727:GSE84727_unmethylatedIntensities.csv.gz
# Post-mortem brain (GSE89702 cerebellum, GSE89703 hippocampus, GSE89705/GSE89706 striatum)
./download_geo.sh raw GSE89702:matrix GSE89703:matrix GSE89705:matrix GSE89706:matrix \
  GSE89702:GSE89702_detectionP_CER_DBCBB.csv.gz GSE89702:GSE89702_methylatedintensities_CER_DBCBB.csv.gz GSE89702:GSE89702_unmethylatedintensities_CER_DBCBB.csv.gz \
  GSE89703:GSE89703_detectionP_HC_LNDBB.csv.gz GSE89703:GSE89703_methylatedintensities_HC_LNDBB.csv.gz GSE89703:GSE89703_unmethylatedintensities_HC_LNDBB.csv.gz \
  GSE89705:GSE89705_detectionP_STR_DBCBB.csv.gz GSE89705:GSE89705_methylatedintensities_STR_DBCBB.csv.gz GSE89705:GSE89705_unmethylatedintensities_STR_DBCBB.csv.gz \
  GSE89706:GSE89706_detectionP_STR_LNDBB.csv.gz GSE89706:GSE89706_methylatedintensities_STR_LNDBB.csv.gz GSE89706:GSE89706_unmethylatedintensities_STR_LNDBB.csv.gz
# Link brain GSM accessions to array IDs (writes sample_sheets/GSE8970x_linked.tsv)
for spec in GSE89702:CER_DBCBB GSE89703:HC_LNDBB GSE89705:STR_DBCBB GSE89706:STR_LNDBB; do
  g=${spec%%:*}; s=${spec#*:}
  Rscript ../inst/scripts/link_geo_samples.R raw/${g}_series_matrix.txt.gz \
    raw/${g}_methylatedintensities_$s.csv.gz raw/${g}_unmethylatedintensities_$s.csv.gz sample_sheets/${g}_linked.tsv
done
# GSE74193 DLPFC IDAT subset (60 samples listed in sample_sheets/GSE74193_subset.tsv; URLs in gse74193_idat_urls.txt)
./download_geo.sh raw GSE74193:matrix
mkdir -p idat/GSE74193 && (cd idat/GSE74193 && xargs -n1 curl -sSLO --retry 5 < ../../gse74193_idat_urls.txt)
# GSE166255 whole-blood RRBS subset (24 samples in sample_sheets/GSE166255_subset.tsv)
./download_geo.sh raw GSE166255:matrix
mkdir -p rrbs/GSE166255 && (cd rrbs/GSE166255 && xargs -n1 curl -sSLO --retry 5 < ../../gse166255_urls.txt)
# Reactome gene sets (offline pathway analysis)
mkdir -p genesets && curl -sL -o genesets/ReactomePathways.gmt.zip https://reactome.org/download/current/ReactomePathways.gmt.zip && (cd genesets && unzip -o ReactomePathways.gmt.zip)
```

Subsets were drawn reproducibly (Python `random.seed(42)` for GSE74193: 30 adult [age >= 18,
dropsample = FALSE] Schizo + 30 Control; `random.seed(7)` for GSE166255: 12 control + 12
neuropathic chronic pain).
