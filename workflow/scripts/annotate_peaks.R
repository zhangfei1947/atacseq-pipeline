suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicRanges)
  library(rtracklayer)
})

dir.create(dirname(snakemake@output[["genes"]]), recursive=TRUE, showWarnings=FALSE)
sink(snakemake@log[[1]], split=TRUE)
bed <- read.delim(snakemake@input[["peaks"]], header=FALSE, stringsAsFactors=FALSE)
colnames(bed)[1:6] <- c("chrom", "start", "end", "peak_id", "score", "strand")
gr <- GRanges(bed$chrom, IRanges(bed$start + 1L, bed$end), peak_id=bed$peak_id, score=bed$score)

# Bioconductor 3.22 moved makeTxDbFromGFF() out of GenomicFeatures into the
# separate txdbmaker package.  The pinned ChIPseeker OCI image does not ship
# txdbmaker, so build the peak-to-gene table directly from the supplied GTF.
# This keeps the annotation tied to the exact reference release and avoids a
# hidden dependency on an organism-specific TxDb package.
gtf <- import(snakemake@input[["gtf"]])
genes <- gtf[gtf$type == "gene"]
if (length(genes) == 0L) {
  stop("No gene features were found in the supplied GTF")
}

gene_id <- if ("gene_id" %in% colnames(mcols(genes))) as.character(mcols(genes)$gene_id) else rep(NA_character_, length(genes))
gene_name <- if ("gene_name" %in% colnames(mcols(genes))) as.character(mcols(genes)$gene_name) else gene_id
gene_tss <- ifelse(as.character(strand(genes)) == "-", end(genes), start(genes))
tss <- GRanges(seqnames(genes), IRanges(gene_tss, gene_tss), strand=strand(genes))

nearest <- distanceToNearest(gr, tss, ignore.strand=TRUE)
nearest_idx <- rep(NA_integer_, length(gr))
nearest_idx[queryHits(nearest)] <- subjectHits(nearest)

peak_center <- floor((start(gr) + end(gr)) / 2)
distance_to_tss <- rep(NA_integer_, length(gr))
valid <- !is.na(nearest_idx)
distance_to_tss[valid] <- peak_center[valid] - gene_tss[nearest_idx[valid]]
minus <- valid & as.character(strand(genes))[nearest_idx] == "-"
distance_to_tss[minus] <- -distance_to_tss[minus]

gene_hits <- findOverlaps(gr, genes, ignore.strand=TRUE)
overlaps_gene <- rep(FALSE, length(gr))
overlaps_gene[unique(queryHits(gene_hits))] <- TRUE
annotation <- ifelse(
  !is.na(distance_to_tss) & distance_to_tss >= -1000L & distance_to_tss <= 100L,
  "Promoter (-1 kb, +100 bp)",
  ifelse(overlaps_gene, "Gene body", "Distal/intergenic")
)

out <- data.frame(
  seqnames=as.character(seqnames(gr)),
  start=start(gr),
  end=end(gr),
  width=width(gr),
  strand=as.character(strand(gr)),
  peak_id=mcols(gr)$peak_id,
  score=mcols(gr)$score,
  annotation=annotation,
  geneChr=ifelse(valid, as.character(seqnames(genes))[nearest_idx], NA_character_),
  geneStart=ifelse(valid, start(genes)[nearest_idx], NA_integer_),
  geneEnd=ifelse(valid, end(genes)[nearest_idx], NA_integer_),
  geneStrand=ifelse(valid, as.character(strand(genes))[nearest_idx], NA_character_),
  geneId=ifelse(valid, gene_id[nearest_idx], NA_character_),
  geneName=ifelse(valid, gene_name[nearest_idx], NA_character_),
  distanceToTSS=distance_to_tss,
  annotation_method="FlyBase GTF nearest TSS; promoter -1000/+100 bp",
  stringsAsFactors=FALSE
)
write.table(out, snakemake@output[["genes"]], sep="\t", quote=FALSE, row.names=FALSE)

if (isTRUE(snakemake@params[["has_redfly"]])) {
  redfly <- import(snakemake@input[["redfly"]][[1]])
  hits <- findOverlaps(gr, redfly, ignore.strand=TRUE)
  evidence <- data.frame(
    peak_id=mcols(gr)$peak_id[queryHits(hits)],
    redfly_seqnames=as.character(seqnames(redfly))[subjectHits(hits)],
    redfly_start=start(redfly)[subjectHits(hits)] - 1L,
    redfly_end=end(redfly)[subjectHits(hits)],
    redfly_name=if ("name" %in% colnames(mcols(redfly))) mcols(redfly)$name[subjectHits(hits)] else "."
  )
  write.table(evidence, snakemake@output[["redfly"]], sep="\t", quote=FALSE, row.names=FALSE)
} else {
  write.table(data.frame(peak_id=character(), redfly_name=character()), snakemake@output[["redfly"]], sep="\t", quote=FALSE, row.names=FALSE)
}
sink()
