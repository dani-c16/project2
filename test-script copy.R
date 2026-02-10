
library(GenomicRanges) # Ranged data library
library(bsseq) # DNA methylation library
library(BSgenome.Hsapiens.UCSC.hg38) # Genome sequence library
library(TxDb.Hsapiens.UCSC.hg38.knownGene) # Gene annotation library
library(Biostrings)
BiocManager::install("dplyr")
library(dplyr)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("qs")
## Get the promoter ranges
 my_promoters <- GenomicRanges::promoters(x = TxDb.Hsapiens.UCSC.hg38.knownGene)
 seqlengths(BSgenome.Hsapiens.UCSC.hg38)
 seqlevels(my_promoters
           )
 geno <- seqlengths(BSgenome.Hsapiens.UCSC.hg38) 
 str(geno)
 geno1 <- seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene)
wholechromo<- GRanges(seqnames = names(geno)
        , ranges = IRanges(start = 1, end = geno)) #gives ranges for whole genome
geno["chr22"]
iran
BSseq()
calc_mC_cov_weighted(bs_obj, gr = wholechromo) #whole chromosome coverage etc.

calc_mC_cov_weighted <- function(RDS_path, gr){
  
  #   message(paste0("Reading ", RDS_path, " ", Sys.time()))
  bs_obj <- readRDS("/Users/dani_c/Downloads/merged_control_CpG-bsseq.Rds")
  
  message("Calculating coverage...")
  gr$Cov <- bsseq::getCoverage(BSseq = bs_obj, regions = gr, type = "Cov",
                               what = "perRegionTotal")
  
  message("Calculating M...")
  gr$M <- bsseq::getCoverage(BSseq = bs_obj, regions = gr, type = "M",
                             what = "perRegionTotal")
  
  message("Calculating methylation percentage...")
  gr$pc <- gr$M / gr$Cov
  
  return(gr)
}

 seqlengths(pro["chr1"])
 pro[""]
 
 seqlengths(geno)
myP <- GenomicRanges::promoters(BSgenome.Hsapiens.UCSC.hg38)
 autosomes <- paste0("chr", c(1:22))
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
 seqnames(my_promoters) %in% autosomes 
 seqnames(pro) %in% autosomes
str(pro
    )
seqlengths(pro)
 pro <- my_promoters[seqnames(my_promoters) %in% autosomes ]
project_data<-readRDS("/Users/dani_c/Downloads/merged_control_CpG-bsseq.Rds")#sample files
my_promoter_sequence <- Biostrings::getSeq(x = BSgenome.Hsapiens.UCSC.hg38,
                                           my_promoters[1:711])
my_promoter_sequence2 <- Biostrings::getSeq(x = BSgenome.Hsapiens.UCSC.hg38,
                                           pro[1:711])
gr <- pro
calc_mC_cov_weighted <- function(RDS_path, gr){
  
  #   message(paste0("Reading ", RDS_path, " ", Sys.time()))
  bs_obj <- readRDS("/Users/dani_c/Downloads/merged_control_CpG-bsseq.Rds")
  
  message("Calculating coverage...")
  gr$Cov <- bsseq::getCoverage(BSseq = bs_obj, regions = gr, type = "Cov",
                               what = "perRegionTotal")
  
  message("Calculating M...")
  gr$M <- bsseq::getCoverage(BSseq = bs_obj, regions = gr, type = "M",
                             what = "perRegionTotal")
  
  message("Calculating methylation percentage...")
  gr$pc <- gr$M / gr$Cov
  
  return(gr)
}

seqnames(my_promoter_sequence2) %in% autosomes 
prom <- my_promoters[seqnames(my_promoter_sequence2) %in% autosomes ]
promo_mean <- calc_mC_cov_weighted(bs_obj, gr = pro)#CG_density <- 
nucleo <-  Biostrings::dinucleotideFrequency(my_promoter_sequence)[ ,"CG"]
 CG <- (nucleo/2200)*100

dinucleotideFrequency(my_promoter_sequence,step = )
Biostrings::dinucleotideFrequency(myBiostrings::dinucleotideFrequency(myBiostrings::dinucleotideFrequency(my_promoter_sequence2)
getCoverage(TxDb.Hsapiens.UCSC.hg38.knownGene, type = M)
merge
table(promo_mean:seqnames, promo_mean:pc)
table(promo_mean)
table(gr)
coverage_values <- getCoverage(bs_obj)
chromosomes <- seqnames(transcripts(txdb)) %in% autosomes
merge
pro2 <- txdb[chromosomes ]
data_frame <- data.frame(
  Chromosome = as.character(gr),
  IRanges = as.character(gr),
  Coverage = calc_mC_cov_weighted( gr = pro)
)
chromo1 <- seqnames(transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)) %in% autosomes
getP
library(dplyr)
data <- promo_mean %>%
  select(column1, column2, column6,column7, column8) 
  select
data.frame2GRanges(data_frame, ignoreStrand = TRUE)
data3 <- data_frame %>%
  select(Chromosome, Coverage.merged, Coverage.merged.1, Coverage.merged.2)
promo_mean[chromosomes]
data4 <- data.frame(promo_mean[chromosomes],
                      Chromosomes = as.character(pro),
                    IRanges = as.character(pro),
                    Coverage = promo_mean) 
loci
  select(seqnames)
ah <- data.seqnames()
ah <- data.frame(promo_mean, 
  Chromosome = as.character(pro),
  IRanges = as.character(promo_mean),
  coverage = promo_mean
)
gr <- GRanges(
  seqnames = Rle(ah$Chromosome),                   # Using Chromosome as seqnames
  ranges = IRanges(start = 1, end = promo_mean)  # Assuming start and end based on promo_mean
)
dat <- ah %>% dplyr::select(seqnames, start, end, merged, merged.1, merged.2)
agg_val <- dat %>%
  group_by(seqnames) %>%
  

print(dat)
summary(dat)
summary(ah)
bs_obj <- readRDS("/Users/dani_c/Downloads/merged_control_CpG-bsseq.Rds")
subset_gr <- gr[seqnames(gr) %in% autosomes]
subset_bsseq <- bs_obj[overlapsAny(gr, subset_gr)]
print(subset_bsseq)
over
gr[seqnames(gr)=="chr1"]
gr2 <- promoters(txdb)
gr2[seqnames(gr2)== "chr1"]
overlapsAny(gr, bs_obj, type=c("any"))
geno1[seqnames(geno1) == "chr1"]
