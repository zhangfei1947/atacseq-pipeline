suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
})

dir.create(dirname(snakemake@output[["genes"]]), recursive=TRUE, showWarnings=FALSE)
sink(snakemake@log[[1]], split=TRUE)
bed <- read.delim(snakemake@input[["peaks"]], header=FALSE, stringsAsFactors=FALSE)
colnames(bed)[1:6] <- c("chrom", "start", "end", "peak_id", "score", "strand")
gr <- GRanges(bed$chrom, IRanges(bed$start + 1L, bed$end), peak_id=bed$peak_id, score=bed$score)
txdb <- makeTxDbFromGFF(snakemake@input[["gtf"]], format="gtf")
anno <- annotatePeak(gr, TxDb=txdb, tssRegion=c(-1000, 100), verbose=FALSE)
out <- as.data.frame(anno)
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
