install.packages("scmeth")
BiocManager::install("scmeth")
library(scmeth)
library(GenomicRanges) # Ranged data library
library(bsseq) # DNA methylation library
help(scmeth)
BiocManager::install("methrix")
library(methrix)
"scmeth"
bs_obj <- readRDS("untitled folder/merged_control_CpG-bsseq.Rds")
bs_seq_obj <- methrix2bsseq(m = bs_obj)
dim(bs_obj)
str(bs_obj)



chromosome_numbers <- paste0("chr",1:22)
chromosomeCoverage(bs_obj[seqnames(values), c(chr1:chr22)])
bs_obj2 <- bs_obj[seqnames(bs_obj) %in% autosomes]
chromosomeCoverage(gr_all2)
chromosomeCoverage(bs_obj2)
data("BS.chr22")
chromosomeCoverage(BS.chr22)

directory <- system.file("extdata/bismark_data",package='scmeth')
bs <- HDF5Array::loadHDF5SummarizedExperiment(directory)
bs2 <- chromosomeCoverage(bs)
table(bs2)

