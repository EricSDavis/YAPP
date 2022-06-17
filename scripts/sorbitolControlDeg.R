## RNA-seq Differential Gene Expression Analysis for hg38

## Load required libraries
library(data.table)
library(tximeta)
library(DESeq2)
library(dplyr)
library(plyranges)
library(liftOver)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

## Import with tximeta ---------------------------------------------------------

## Read in sample sheet as coldata
coldata <- 
  file.path("rawdata", "rna", "YAPP_HEK_wt_1_RNApipeSamplesheet.txt") |>
  fread() |>
  as.data.frame(stringsAsFactors = TRUE)

## Add quant paths and names
coldata$files <- file.path("rawdata", "rna", "quant", coldata$sn, "quant.sf")
colnames(coldata) <- gsub("sn", "names", colnames(coldata))

## Check that quant paths are correct
file.exists(coldata$files)

## Import data with tximeta & summarize to gene
se <- tximeta(coldata)
gse <- summarizeToGene(se)


## Differential Expression Analysis --------------------------------------------

## Convert to factors (avoids a warning)
colData(gse)[] <- lapply(colData(gse), factor)

## Build DESeq object
dds <- DESeqDataSet(gse, design = ~Bio_Rep + Condition)

## Filter out lowly expressed genes
## (at least 10 counts in at least 2 samples)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep,]

## Fit model
dds <- DESeq(dds)

## Get DEG results as GRanges
de_genes <- 
  results(dds,
          name = "Condition_Sorb_vs_Control",
          format = "GRanges") |>
  plyranges::names_to_column("gene_id")


## Add gene symbols to de_genes and
## remove non-standard chromosomes
de_genes <- 
  inner_join(x = as.data.table(de_genes),
             y = as.data.table(rowData(gse)) |>
               dplyr::select(c("gene_id", "symbol", "tx_ids")),
             by = "gene_id") |>
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) |>
  keepStandardChromosomes(pruning.mode = "coarse")

## Add seqinfo to object
txdb <- 
  TxDb.Hsapiens.UCSC.hg38.knownGene |>
  keepStandardChromosomes()

seqlevels(de_genes) <- seqlevels(txdb)
seqinfo(de_genes) <- seqinfo(txdb)


## Lift over to hg19 -----------------------------------------------------------

## Create new object and change seqlevelstyles
de_genes_hg19 <- de_genes

## Import chain
ch <-
  system.file("extdata", "hg38ToHg19.over.chain",
              package = "liftOver") |>
  import.chain()

## Lift over and set genome
de_genes_hg19 <- 
  liftOver(de_genes_hg19, ch) |>
  unlist()

## Add seqinfo to object
txdb <- 
  TxDb.Hsapiens.UCSC.hg19.knownGene |>
  keepStandardChromosomes()

seqlevels(de_genes_hg19) <- seqlevels(txdb)
seqinfo(de_genes_hg19) <- seqinfo(txdb)


## Save objects ----------------------------------------------------------------

saveRDS(de_genes, file = file.path("data", "sorbitolControlDegHg38.rds"))
saveRDS(de_genes_hg19, file = file.path("data", "sorbitolControlDegHg19.rds"))
