install.packages(GenomicRanges)
BiocManager::install("GenomicRanges")
BiocManager::install("Biostrings")
BiocManager::install("scmeth")
BiocManager::install("bsseq")

library(bsseq)
library(GenomicRanges)
seqlengths(hg38)
autosomes <- paste0("chr", c(1:22))
seqlengths(hg38) %in% autosomes

seqlengths(h38)
seqn
calc_mC_cov_weighted <- function(RDS_path, gr){
  
  #   message(paste0("Reading ", RDS_path, " ", Sys.time()))
  bs_obj <- readRDS("/Users/dani_c/Downloads/merged_control_CpG-bsseq.Rds")
  
  message("Calculating coverage...")
  gr$Cov <- bsseq::getCoverage(BSseq = bs_obj, regions = gr, type = "Cov",
                               what = "perRegionAverage")
  
  message("Calculating M...")
  gr$M <- bsseq::getCoverage(BSseq = bs_obj, regions = gr, type = "M",
                             what = "perRegionAverage")
  
  message("Calculating methylation percentage...")
  gr$pc <- gr$M / gr$Cov
  
  return(gr)
}
print(gr$M)
print(gr$Cov)
wholechromo<- GRanges(seqnames = names(geno)
                      , ranges = IRanges(start = 1, end = geno))
wholechromo1 <- GRanges(seqnames = names(geno)
                      , ranges = IRanges(start = 1, end = geno))
wholechromo[1:12]

calc_mC_cov_unweighted <- function(RDS_path, wholechromo1){
  bs_obj <- readRDS("/Users/dani_c/Downloads/merged_control_CpG-bsseq.Rds")
  
  message("Calculating M...")
  wholechromo1$M <- bsseq::getCoverage(BSseq = bs_obj, regions = wholechromo1, type = "M", what = "perRegionTotal")
  
  message("Calculating methylation coverage...")
  wholechromo1$Cov <- bsseq::getCoverage(BSseq = bs_obj, regions = wholechromo1, type = "Cov", what = "perRegionTotal")
  
  message("Calculating methylation percentage ...")
 wholechromo1$pc <- wholechromo1$M / wholechromo1$Cov
 
 return(wholechromo1)
  }
wholechromo1[1:12]

Mchromo <- getCoverage(bs_obj, regions = wholechromo[1:12], type = "M", what = "perRegionTotal")
covchromo <- getCoverage(bs_obj, regions = wholechromo[1:12], type = "Cov", what = "perRegionTotal")
pcchromo <- Mchromo/covchromo
base::table()
unweightedchromo <- data.frame(
  Chromosome = as.character(wholechromo[1:12]),
  IRanges = as.character(wholechromo[1:12]),
  M = Mchromo,
  Coverage = covchromo,
  PC = pcchromo
)


wholechromo<- GRanges(seqnames = names(geno)
        , ranges = IRanges(start = 1, end = geno["chr22"]))
getCoverage(BSseq = bs_obj, regions = wholechromo[1:12], type = "M", what = "perRegionTotal")
