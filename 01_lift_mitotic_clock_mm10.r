#Libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Mmusculus.UCSC.mm10)

library(liftOver)
#folder
data.dir <- "/home/epicwl/c010-datasets/Internal/mitotic_clock/data/"

#load cpgs
cpg_human <- read.table(file.path(data.dir, "cpgs_human.csv"), sep=",", header=TRUE)
cpg_human<- cpg_human[, c("ID", "Chromosome", "Start","End", "cgid")]

cpg_human_gr <- GRanges(
  seqnames = cpg_human$Chromosome,
  ranges = IRanges(start = cpg_human$Start,
                   end = cpg_human$End
  )
)
mcols(cpg_human_gr)<- cpg_human$cgid

#annotate CpGs
txdb_hs <- TxDb.Hsapiens.UCSC.hg19.knownGene
anno_cpg_human <- annotatePeak(peak= cpg_human_gr, tssRegion=c(-1500, 1500),
                         TxDb=txdb_hs, annoDb="org.Hs.eg.db")
anno_cpg_human

#get liftover
ch <- import.chain(file.path(data.dir, "hg19ToMm10.over.chain"))
seqlevelsStyle(cpg_human_gr) = "UCSC"
cpg_mouse_gr<-liftOver(cpg_human_gr, ch)

#remove the ones which were split up or are 0 
table(unlist(lapply(cpg_mouse_gr, function(x)length(x))))
idx <- unlist(lapply(cpg_mouse_gr, function(x)length(x))) == 1
#subset list
cpg_mouse_gr_sub <- cpg_mouse_gr[idx]
#unlist them
cpg_mouse_gr_sub<- unlist(cpg_mouse_gr_sub)

#annotate new mouse cpgs
txdb_mm <- TxDb.Mmusculus.UCSC.mm10.knownGene
anno_cpg_mouse <- annotatePeak(peak= cpg_mouse_gr_sub, tssRegion=c(-1500, 1500),
                         TxDb=txdb_mm, annoDb="org.Mm.eg.db")


#find out if they are cpg
mm10Genome <- BSgenome.Mmusculus.UCSC.mm10
mouse_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10,cpg_mouse_gr_sub )
table(mouse_seq)

#inex the ones which are cpgs 
idx_cg <- mouse_seq=="CG"
#subset the ones which 
cpg_mouse_gr_sub_cg <- cpg_mouse_gr_sub[idx_cg,]
#save them as a bed file
export.bed(cpg_mouse_gr_sub_cg, file.path(data.dir, "cpgs_lifted_mm10_sub.bed"))

#get promoter regions of genes
mm10_promoters <- promoters(genes(txdb_mm), upstream= 1500, downstream=1500)
anno_cpg_mouse_df <- as.data.frame(anno_cpg_mouse)
anno_cpg_mouse_df$cpgid <- cpg_mouse_gr_sub$X
anno_cpg_mouse_df_unique <- anno_cpg_mouse_df[duplicated(anno_cpg_mouse_df$geneID),]
relevant_promoters_mm10 <- mm10_promoters[mm10_promoters$gene_id %in% as.data.frame(anno_cpg_mouse)$geneId, ]
relevant_promoters_mm10$symbols <- unlist(mapIds(org.Mm.eg.db,relevant_promoters_mm10$name, "SYMBOL", "ENTREZID", multiVals = "first"))
export.bed(relevant_promoters_mm10, file.path(data.dir, "promoters_relCpgs_mm10.bed"))

´
