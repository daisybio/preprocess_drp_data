library(SummarizedExperiment)
library(tximport)
library(data.table)

args <- c("/path/to/ensembl107_GRCh38/Homo_sapiens.GRCh38.107.gtf",
          "/path/to/nf_core/star_salmon",
          "ccle_salmon")

path <- args[2]
prefix <- args[3]
run_table <- fread("/path/to/SraRunTable_CCLE.txt")
mapping <- run_table[, c("Run", "Cell_line")]
mapping$Cell_line[is.na(mapping$Cell_line)] <- 'T.T'
mapping$Cell_line[mapping$Cell_line == 'HARA'] <- 'HARA [Human squamous cell lung carcinoma]'
mapping$Cell_line[mapping$Cell_line == 'ML-1'] <- 'ML-1 [Human thyroid carcinoma]'
mapping$Cell_line[mapping$Cell_line == 'KS-1'] <- 'KS-1 [Human glioblastoma]'
mapping$Cell_line[mapping$Cell_line == 'HH'] <- 'HH [Human lymphoma]'
mapping$Cell_line[mapping$Cell_line == 'RCM-1'] <- 'RCM-1 [Human rectal adenocarcinoma]'
mapping$Cell_line[mapping$Cell_line == 'G-292\\, clone A141B1'] <- 'G-292 clone A141B1'
mapping$Cell_line[mapping$Cell_line == 'CJM'] <- 'CJM [Human melanoma]'
mapping$Cell_line[mapping$Cell_line == 'KD'] <- 'KD [Human abdomen rhabdoid tumor]'

cellosaurus_ids <- fread('/path/to/cell_line_mapping.csv', header=T)

# reformat SY:
cellosaurus_ids[, SY := toupper(gsub("'", "", gsub(",\\s*", ",", gsub("(\\[|\\])", "", SY))))]

get_ac_value <- function(query_id){
  query_id <- toupper(query_id)
  query_id <- unlist(strsplit(query_id,';'))[1]
  cid <- cellosaurus_ids[query_id == toupper(ID), AC]
  if(identical(cid, character(0))){
    ids_list <- strsplit(cellosaurus_ids$SY, ',\\s*')
    matching_rows <- which(sapply(ids_list, function(ids) query_id %in% ids))
    if(length(matching_rows) > 0){
      cid <- cellosaurus_ids$AC[matching_rows[[1]]]
    }else{
      cid <- ""
      print(query_id)
    }
  }else if(length(cid) == 2){
    cid <- cellosaurus_ids[query_id == toupper(ID), AC][[1]]
  }
  return(cid)
}

cid_list <- sapply(mapping$Cell_line, get_ac_value)
mapping[, Cellosaurus_ID := cid_list[Cell_line]]
fwrite(mapping, '/path/to/mapping_Sra_names_cellosaurus.csv')

tx2gene <- "/path/to/nf_core/star_salmon/salmon_tx2gene.tsv"
info <- file.info(tx2gene)
if (info$size == 0) {
  tx2gene <- NULL
} else {
  rowdata <- read.csv(tx2gene, sep="\t", header = FALSE)
  colnames(rowdata) <- c("tx", "gene_id", "gene_name")
  tx2gene <- rowdata[,1:2]
}

fns <- list.files(path, pattern = "quant.sf", recursive = T, full.names = T)
names <- basename(dirname(fns))
names(fns) <- names

coldata <- data.frame(files = fns, names = names)
coldata <- merge(coldata, mapping, by.x = "names", by.y = "Run")

txi <- tximport(fns, type = "salmon", txOut = TRUE)
rownames(coldata) <- coldata[["names"]]
extra <- setdiff(rownames(txi[[1]]),  as.character(rowdata[["tx"]]))
if (length(extra) > 0) {
  rowdata = rbind(rowdata, data.frame(tx=extra, gene_id=extra, gene_name=extra))
}
rowdata <- rowdata[match(rownames(txi[[1]]), as.character(rowdata[["tx"]])),]
rownames(rowdata) <- rowdata[["tx"]]
se <- SummarizedExperiment(assays = list(counts = txi[["counts"]], abundance = txi[["abundance"]], length = txi[["length"]]),
                          colData = DataFrame(coldata),
                          rowData = rowdata)
if (!is.null(tx2gene)) {
  gi <- summarizeToGene(txi, tx2gene = tx2gene)
  gi.ls <- summarizeToGene(txi, tx2gene = tx2gene,countsFromAbundance="lengthScaledTPM")
  gi.s <- summarizeToGene(txi, tx2gene = tx2gene,countsFromAbundance="scaledTPM")
  growdata <- unique(rowdata[,2:3])
  growdata <- growdata[match(rownames(gi[[1]]), growdata[["gene_id"]]),]
  rownames(growdata) <- growdata[["tx"]]
  gse <- SummarizedExperiment(assays = list(counts = gi[["counts"]], abundance = gi[["abundance"]], length = gi[["length"]]),
                             colData = DataFrame(coldata),
                             rowData = growdata)
  gse.ls <- SummarizedExperiment(assays = list(counts = gi.ls[["counts"]], abundance = gi.ls[["abundance"]], length = gi.ls[["length"]]),
                                colData = DataFrame(coldata),
                                rowData = growdata)
  gse.s <- SummarizedExperiment(assays = list(counts = gi.s[["counts"]], abundance = gi.s[["abundance"]], length = gi.s[["length"]]),
                               colData = DataFrame(coldata),
                               rowData = growdata)
}

build_table <- function(se.obj, slot) {
  tmp <- cbind(rowData(se.obj)[,1:2], assays(se.obj)[[slot]])
  match_colnames <- data.table(names = colnames(tmp[, -c(1,2)]))
  match_colnames <- merge(match_colnames, colData(se.obj))
  colnames(tmp) <- c(names(rowData(se.obj)[,1:2]), match_colnames$Cellosaurus_ID)
  return(tmp)
}

if(exists("gse")){
  write.table(build_table(gse, "abundance"), paste(c(prefix, "gene_tpm_cvcl.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
  write.table(build_table(gse, "counts"), paste(c(prefix, "gene_counts_cvcl.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
  write.table(build_table(gse, "length"), paste(c(prefix, "gene_length_cvcl.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
  #write.table(build_table(gse.ls, "abundance"), paste(c(prefix, "gene_tpm_length_scaled_cvcl.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
  #write.table(build_table(gse.ls, "counts"), paste(c(prefix, "gene_counts_length_scaled_cvcl.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
  #write.table(build_table(gse.s, "abundance"), paste(c(prefix, "gene_tpm_scaled_cvcl.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
  write.table(build_table(gse.s, "counts"), paste(c(prefix, "gene_counts_scaled_cvcl.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
}

write.table(build_table(se,"abundance"), paste(c(prefix, "transcript_tpm_cvcl.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(build_table(se, "counts"), paste(c(prefix, "transcript_counts_cvcl.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
write.table(build_table(se, "length"), paste(c(prefix, "transcript_length_cvcl.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)

# Print sessioninfo to standard out
citation("tximeta")
sessionInfo()