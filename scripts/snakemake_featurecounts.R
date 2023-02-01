library(Rsubread)
library(tibble)
library(readr)
library(dplyr)

reference_df <- data.frame(reference = c("NC_000964.3", "NZ_CP014269.1", "NC_000866.4", "NC_011421.1"),
                           name = c("B. subtilis", "E. coli", "T4 phage", "SPO1 phage"), 
                           class = c("Bacteria", "Bacteria", "Phage", "Phage"))

fc <- featureCounts(files = snakemake@input[['bam']],
                    isGTFAnnotationFile = T,
                    annot.ext = snakemake@input[['gtf']],
                    isPairedEnd = T,
                    GTF.featureType = "gene")

fc_counts_tmp <- fc$counts %>%
  as.data.frame() %>%
  tibble::rownames_to_column("GeneID")

fc_counts <- left_join(reference_df, as.data.frame(fc$annotation), by = c("reference" = "Chr"), multiple  = "all") %>%
  left_join(fc_counts_tmp, by = "GeneID")

write_tsv(fc$stat, snakemake@output[['stat']])
write_tsv(fc_counts, snakemake@output[['counts']])