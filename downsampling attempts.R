install.packages("DropletUtils")
BiocManager::install("DropletUtils")
library(DropletUtils)
library(bsseq)
library(GenomicRanges)
library(Biostrings) 
library(scmeth)
BiocManager::install("rhdf5")
library(rhdf5)
h5file <- h5createFile(h5mergedreads)
bs_obj <- readRDS('/Users/dani_c/Downloads/merged_control_CpG-bsseq copy.Rds')
#bs_obj2 <- seqnames(bs_obj) %in% autosomes

# Extract methylation levels
methylation_levels <- getMeth(bs_obj, type = "raw")
coverage_levels <- getCoverage(bs_obj)
chromosomes_rle <- seqnames(bs_obj)
positions <- start(bs_obj)
chromosomes <- as.vector(as.character(chromosomes_rle)) 
# You may want to convert to a data frame for easier saving
data_to_save <- data.frame(
  chromosome = as.character(chromosomes),
  position = positions,
  methylation = as.vector(methylation_levels),     # Convert to a vector if needed
  coverage = as.vector(coverage_levels)             # Convert to a vector if needed
)


h5_file_path <- "output_bsseq_file2.h5"
h5createFile(h5_file_path) 


h5write(methylation_levels,file = "h5_file_path", name = "methylation_levels")
h5write(coverage_levels, file = "h5_file_path", name = "coverage_levels")
h5write(data_to_save, h5_file_path, "data_frame")
h5write(methylation_levels, h5_file_path, "methylation_levels")

h5write(positions,h5_file_path, "positions")
h5write(chromosomes, file = h5_file_path, "Chromosomes")
methyl_verified <- h5read(h5_file_path, "methylation_levels")

wholc<-as.data.frame(wholechromo2[1:22])


as.matrix.data.frame(wholc)
graph[,"chr"1:22)]
wholechromo2[1:22]

downsampleReads(h5_file_path,prop = 0.5,barcode.length = NULL, bycol = FALSE)
downsampleReads(gr,prop = 0.5,barcode.length = NULL, bycol = FALSE)
ds_cov_df <- downsampleMatrix(coverage_df, prop = 0.5, bycol = TRUE)
downsampleMatrix(my_TF, prop = 0.5, bycol =  TRUE)


seqlengths(ds_cov_df)

h5ls("/Users/dani_c/Documents/GitHub/project/output_bsseq_file2.h5")

downsample(bs_obj, dsRates = 0.5, subSample = 1e+6, offset = 0)
directory <- system.file("extdata/bismark_data", package='scmeth')
bs <- HDF5Array::loadHDF5SummarizedExperiment(directory)
scmeth::downsample(bs)
barc
