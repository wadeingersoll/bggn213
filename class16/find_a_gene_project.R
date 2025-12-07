library(bio3d)

msa <- read.fasta("/Users/wadeingersoll/Downloads/aligned_sequences.fst")

matrix <- seqidentity(msa$ali)

heatmap(matrix)

png("heatmap.png", width = 1000, height = 1000, res = 150)
heatmap(matrix,
        scale = "none",
        margins = c(12, 12),
        cexRow = 1.1, cexCol = 1.1)
dev.off()

## Finding sequence with highest identity to all others

msa <- read.fasta("/Users/wadeingersoll/Downloads/aligned_sequences.fst")
matrix <- seqidentity(msa$ali)

diag(matrix) <- NA  # ignore self-identity (always 1)
row_avg <- rowMeans(matrix, na.rm = TRUE)
row_max <- apply(matrix, 1, max, na.rm = TRUE)

idx_avg <- which.max(row_avg)   # by average similarity
idx_max <- which.max(row_max)   # by single highest similarity (if you prefer)

# See which sequence it is
msa$id[idx_avg]   # the ID (name) of your chosen sequence

query <- paste(msa$ali[idx_avg, ], collapse="")
query <- gsub("\\t", "", query)

blast <- blast.pdb(query)
annotated_hit_table <- pdb.annotate(blast$hit.tbl$subjectids)

hit_table_df <- data.frame(ID=annotated_hit_table$structureId,
                           Technique=annotated_hit_table$experimentalTechnique,
                           Resolution=annotated_hit_table$resolution,
                           Source=annotated_hit_table$source,
                           Evalue=blast$hit.tbl$evalue,
                           Identity=blast$hit.tbl$identity)

