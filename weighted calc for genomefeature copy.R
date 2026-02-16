install.packages("writexl")
library(writexl)
library(GenomicRanges)
my_promoters
chromo <- wholechromo[1:22]
my_terminators <- GenomicRanges::terminators(x = TxDb.Hsapiens.UCSC.hg38.knownGene, c(1:22))
#separate out the 22 chromosomes
my_terminator <- my_terminators[seqnames(my_terminators) %in% autosomes]
my_promoters
#for chromosome 1 of hg38
my_enhancers <- read.csv('/Users/dani_c/Desktop/enhancer ranges for hg38 2.csv')
#make GRanges
gr_enhance <- makeGRangesFromDataFrame(my_enhancers, 
                               keep.extra.columns = TRUE, 
                               seqnames.field = "seqnames", 
                               start.field = "start", 
                               end.field = "end", 
                             )
#weighted calc for only chromosome 1 enhancers of 100X sample
enh <- calc_mC_cov_weighted(bs_obj, gr_enhance) |> as.data.frame(
  
)

enh[1:100,]

#terminators cov adn M for 22 chromosomes
calc_mC_cov_weighted(bs_obj,my_terminator)
#separating out chromosome 1
chr1 <- paste0("chr", c(1))
my_term <- my_terminators[seqnames(my_terminators) %in% chr1]
terminators <- calc_mC_cov_weighted(bs_obj, my_term)
summary(terminators)
print(terminators)
#now as an excel file
graph_weighted <- as.data.frame(terminators)
t_uw <- calc_mC_cov_unweighted(bs_obj, my_term)
uw_term <- as.data.frame(t_uw)
write_xlsx(graph, "output.xlsx")
cbind
df3 <- data.frame(graph_weighted$merged.2, uw_term$merged.2)

#Unit Test -Esque

M <- matrix(c(10,10,5),3,1)
Cov <- matrix(c(50,100,10),3,1)
BS1 <- BSseq(chr = c("chr1","chr1", "chr1"), pos = c(1,2,3), M = M, Cov = Cov, sampleNames = c("A"))
GR1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1,end = 3))
calc_mC_cov_unweighted(BS1, GR1 )
calc_mC_cov_weighted(BS1, GR1)
c <- getCoverage(BS1, regions = GR1, type = "Cov", what = "perRegionAverage")
d <- getCoverage(BS1, regions = GR1, type = "M", what = "perRegionAverage")


a <- getCoverage(BS1, regions = GR1, type = "Cov", what = "perRegionTotal")
b <- getCoverage(BS1, regions = GR1, type = "M", what = "perRegionTotal")
b/a
d/c



#unweighted
getMeth(BS1, regions = GR1, type = "raw", what = "perRegion")



uw_true <- getMeth(bs_obj, regions = pro, type = "raw", what = "perRegion")
cov_w_pro <- calc_mC_cov_weighted(bs_obj, pro)
pc_uw <- as.data.frame(uw_true)
df_cov_w_pro <- as.data.frame(cov_w_pro)
df4 <- data.frame(df_cov_w_pro, pc_uw
                  )
df4$delta <- abs(df4$merged.2 - df4$merged.3)
plot(df4$merged, df4$delta, xlim = c(0,50000))
hist(df4$delta)
hist(df4$merged[df4$merged < 3500])
table()
poissonGoodnessOfFit(BS1, nQuantiles = 3^1)
poissonGoodnessOfFit(bs_obj, nQuantiles = 200^5)
binomialGoodnessOfFit(bs_obj, method = )