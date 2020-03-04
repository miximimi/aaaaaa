


snps.bed.file <- "motifbreakR_SFARI1000_SNVs_hg38_forMPILEUP.bed"
snps.mb <- snps.from.file(snps.bed.file,
                          search.genome = BSgenome.Hsapiens.UCSC.hg38,
                          format = "bed")

pca.enhancer.snps <- sample(snps.mb, 5)


data(hocomoco)
motifs <- hocomoco#sample(hocomoco, 50)

results <- motifbreakR(pca.enhancer.snps,
                       motifs, threshold = 0.85,
                       method = "ic",
                       BPPARAM=BiocParallel::SerialParam())

plotMB(results, "chr5:22065024:A:T", reverseMotif = TRUE, effect = c("strong",
                                                      "weak"))

calculatePvalue(results, background = c(A = 0.25, C = 0.25, G = 0.25, T =
                                          0.25))

### given example result ###
example_result <- example.results[names(example.results) %in% "rs2661839"]
calculatePvalue(example_result)
plotMB(example.results, "rs2661839", reverseMotif = TRUE, effect = c("strong",
  
                                                                                                                                       "weak"))
plotMB(example.results, "rs2661839", effect = "strong")
plotMB(example.results, "rs7012442", effect = "weak")
plotMB(example.results, "rs10486567", effect = "strong") #previous publication cancer

example_result <- example.results[names(example.results) %in% "rs10486567"]
calculatePvalue(example_result)

data(example.results)
example.results
## Not run:
plotMB(example.results, "rs2661839", effect = "strong")
## End(Not run)



# ==============
# try motifbreakr

rs1006140

library("BSgenome.Hsapiens.NCBI.GRCh38")
BiocManager::install("SNPlocs.Hsapiens.dbSNP150.GRCh38")

a <- getSeq(Hsapiens,"chr19")

snps.mb <- snps.from.rsid(rsid = "rs1006140",
                          dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,
                          search.genome = BSgenome.Hsapiens.NCBI.GRCh38)



