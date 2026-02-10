library(GenomicRanges)
seqlengths(geno)
#chromo_data ,- data.frame(
#  seqnames = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
               "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
 # start = c(1, 
            
calc_mC_cov_weighted <- function(RDS_path, wholechromo1){
  
  #   message(paste0("Reading ", RDS_path, " ", Sys.time()))
  bs_obj <- readRDS("/Users/dani_c/Downloads/merged_control_CpG-bsseq.Rds")
  
  message("Calculating coverage...")
  wholechromo1$Cov <- bsseq::getCoverage(BSseq = bs_obj, regions = wholechromo1, type = "Cov",
                               what = "perRegionAverage")
  
  message("Calculating M...")
  wholechromo1$M <- bsseq::getCoverage(BSseq = bs_obj, regions = wholechromo1, type = "M",
                             what = "perRegionAverage")
  
  message("Calculating methylation percentage...")
  wholechromo1$pc <- wholechromo1$M / wholechromo1$Cov
  
  return(wholechromo1)
}