## Sam's function to calculate methylation levels for 1 sample. Function returns a GRanges object
## Function arguments
        # qs_path = path to qs file that has a bsseq object with DNA methylation data
        # gr = Genomic ranges object for ranges to calculate DNA methyaltion levels

## Function to calculate mC levels for one sample.
        ## Returns a granges object
calc_mC_cov_weighted <- function(RDS_path, gr){
        
     #   message(paste0("Reading ", RDS_path, " ", Sys.time()))
        bs_obj <- readRDS(RDS_path)
        
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

message("Calculating M...")
gr$M <- bsseq::getCoverage(BSseq = bs_obj, regions = gr, type = "M",
                           what = "perBase")
