suppressPackageStartupMessages(library(DESeq2))

dir.create(dirname(snakemake@output[["results"]]), recursive=TRUE, showWarnings=FALSE)
sink(snakemake@log[[1]], split=TRUE)

tab <- read.delim(snakemake@input[["counts"]], check.names=FALSE)
meta <- read.delim(snakemake@input[["manifest"]], check.names=FALSE, stringsAsFactors=FALSE)
sample_ids <- meta$sample_id
count_mat <- as.matrix(tab[, sample_ids, drop=FALSE])
rownames(count_mat) <- tab$peak_id
storage.mode(count_mat) <- "integer"

factor_name <- snakemake@params[["factor"]]
numerator <- as.character(snakemake@params[["numerator"]])
denominator <- as.character(snakemake@params[["denominator"]])
meta[[factor_name]] <- relevel(factor(meta[[factor_name]]), ref=denominator)
rownames(meta) <- meta$sample_id
meta <- meta[colnames(count_mat), , drop=FALSE]

design_formula <- as.formula(snakemake@params[["design"]])
mm <- model.matrix(design_formula, data=meta)
if (qr(mm)$rank < ncol(mm)) stop("Design matrix is not full rank")

dds <- DESeqDataSetFromMatrix(countData=count_mat, colData=meta, design=design_formula)
level_n <- table(meta[[factor_name]])
minimum_samples <- min(level_n[c(numerator, denominator)])
keep <- rowSums(counts(dds) >= as.integer(snakemake@params[["min_count"]])) >= minimum_samples
dds <- dds[keep,]
if (nrow(dds) == 0) stop("No peaks remain after independent count prefilter")
dds <- DESeq(dds, parallel=FALSE)

contrast <- c(factor_name, numerator, denominator)
res <- results(dds, contrast=contrast, alpha=as.numeric(snakemake@params[["fdr"]]))
shrunk <- lfcShrink(dds, contrast=contrast, res=res, type="normal")
out <- as.data.frame(res)
out$shrunken_log2FoldChange <- shrunk$log2FoldChange
out$peak_id <- rownames(out)
coords <- tab[match(out$peak_id, tab$peak_id), c("chrom", "start_1based", "end")]
coords$start <- coords$start_1based - 1L
out <- cbind(coords[, c("chrom", "start", "end")], out)
out <- out[, c("peak_id", "chrom", "start", "end", setdiff(colnames(out), c("peak_id", "chrom", "start", "end")))]
write.table(out, snakemake@output[["results"]], sep="\t", quote=FALSE, row.names=FALSE)

fdr <- as.numeric(snakemake@params[["fdr"]])
lfc <- as.numeric(snakemake@params[["lfc"]])
opening <- out[!is.na(out$padj) & out$padj <= fdr & out$log2FoldChange >= lfc,]
closing <- out[!is.na(out$padj) & out$padj <= fdr & out$log2FoldChange <= -lfc,]
write_bed <- function(x, path) {
  if (nrow(x) == 0) { file.create(path); return() }
  bed <- data.frame(x$chrom, x$start, x$end, x$peak_id, x$log2FoldChange, ".")
  write.table(bed, path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
write_bed(opening, snakemake@output[["opening"]])
write_bed(closing, snakemake@output[["closing"]])

norm <- as.data.frame(counts(dds, normalized=TRUE))
norm$peak_id <- rownames(norm)
write.table(norm[, c("peak_id", sample_ids)], snakemake@output[["normalized"]], sep="\t", quote=FALSE, row.names=FALSE)

vsd <- vst(dds, blind=FALSE)
pca <- prcomp(t(assay(vsd)))
pdf(snakemake@output[["pca"]], width=7, height=6)
plot(pca$x[,1], pca$x[,2], pch=19, col=as.integer(meta[[factor_name]]), xlab="PC1", ylab="PC2")
text(pca$x[,1], pca$x[,2], labels=rownames(pca$x), pos=3, cex=0.7)
legend("topright", legend=levels(meta[[factor_name]]), col=seq_along(levels(meta[[factor_name]])), pch=19)
dev.off()

summary_json <- sprintf(
  '{\n  "tested_peaks": %d,\n  "opening_peaks": %d,\n  "closing_peaks": %d,\n  "fdr": %g,\n  "abs_log2fc": %g,\n  "normalization_warning": "Compositional normalization does not establish an absolute genome-wide accessibility shift."\n}\n',
  nrow(out), nrow(opening), nrow(closing), fdr, lfc
)
writeLines(summary_json, snakemake@output[["summary"]])
sink()
