#!/usr/bin/env Rscript
# Build a sample sheet linking GEO GSM accessions to raw-data columns by methylation profile.
# Usage: Rscript link_geo_samples.R <series_matrix.txt.gz> <meth.csv.gz> <unmeth.csv.gz> <out_sample_sheet.tsv>
a <- commandArgs(TRUE)
if (length(a) != 4) stop("Usage: link_geo_samples.R series_matrix meth_file unmeth_file out.tsv")
suppressPackageStartupMessages(library(methena))
lk <- mf_geo_link_by_profile(a[1], a[2], a[3])
ph <- mf_read_geo_pheno(a[1])
ss <- merge(ph, lk, by = "geo_accession")
ss$data_id[!ss$accepted] <- NA
mf_write_tsv(ss, a[4])
if (any(!ss$accepted)) message("WARNING: ", sum(!ss$accepted), " samples could not be linked unambiguously and have data_id = NA")
message("Wrote ", a[4])
