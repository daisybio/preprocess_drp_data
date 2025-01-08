library(data.table)
library(stringr)

mapping <- fread('../../mapping/cellosaurus_01_2024.csv')
cell_line_names <- fread('cell_line_names.csv')
mapping_short <- mapping[, c("AC", "DR")]
# in DR, in the whole string find only the DepMap IDs "ACH-\d{6}" and write it to a new column
mapping_short[, depmap := str_extract(DR, "ACH-\\d{6}")]
mapping_short <- merge(mapping_short, cell_line_names, by.x = "AC", by.y = "cellosaurus_id")
mapping_short <- mapping_short[, -c("DR")]

mutations <- fread('CCLE_mutations.csv')
away <- c("Intron", "Silent", "5'Flank", "3'UTR", "5'UTR")
mutations <- mutations[!Variant_Classification %in% away]
mutations <- mutations[, c("Hugo_Symbol", "DepMap_ID")]
mutations <- merge(mutations, mapping_short, by.x = "DepMap_ID", by.y = "depmap")
mutations <- mutations[, -c("DepMap_ID")]
mutations <- unique(mutations)
mutations$mutated <- "True"
colnames(mutations) <- c("Hugo_Symbol", "cellosaurus_id", "cell_line_name", "mutated")
mutations <- dcast(mutations, cellosaurus_id + cell_line_name ~ Hugo_Symbol, value.var = "mutated")
# fill NA's with "False"
mutations[is.na(mutations)] <- "False"
fwrite(mutations, 'mutations_cellosaurus.csv')
