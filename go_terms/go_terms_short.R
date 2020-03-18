library(miRNAtap)
library(miRNAtap.db)
library(topGO)
library(org.Hs.eg.db)
library(annotate)
library(readxl)

setwd("I:/")
# Evaluate whether vector is empty or not
# https://stackoverflow.com/a/29135732
vector.is.empty <- function(x) return(length(x) == 0)

mir <- read_excel("miRNAs_test.xlsx", col_names = FALSE)
colnames(mir) <- c("miRNAs")
mir <- mir$miRNAs
# Get the most frequent GO terms associated with the miRNAs
# https://bioconductor.org/packages/devel/bioc/vignettes/miRNAtap/inst/doc/miRNAtap.pdf
mirs_data <- data.frame(GO.ID = character(), Term = character(), Annotated = numeric(),
                       Significant = numeric(), Expected = numeric(), KS = character(),
                       stringsAsFactors = FALSE)
for (i in 1:length(mir)){
  print(mir[i])
  predictions <- getPredictedTargets(mir[i], species = "hsa", method = "geom", min_src = 2)
  rankedGenes <- predictions[,"rank_product"]
  selection <- function(x) TRUE
  allGO2genes <- annFUN.org(whichOnto = "BP", feasibleGenes = NULL, mapping = "org.Hs.eg.db", ID = "entrez")
  GOdata <- new("topGOdata", ontology = "BP", allGenes = rankedGenes,
                annot = annFUN.GO2genes, GO2genes = allGO2genes, geneSel = selection, nodeSize=5)
  results.ks <- runTest(GOdata, algorithm = "classic", statistic = "ks")
  allRes <- GenTable(GOdata, KS = results.ks, orderBy = "KS", numChar = 10000)
  mirs_data <- rbind(mirs_data, allRes)
}

rm(i, predictions, rankedGenes, selection, allGO2genes, GOdata, results.ks, allRes)

# Count the GO terms and filter the less frequent ones
mirs_dataTerm <- as.data.frame(mirs_data$Term)
mirs_dataTerm <- as.data.frame(table(mirs_dataTerm[1]))
colnames(mirs_dataTerm)[1] <- c("Term")
mirs_dataTerm <- mirs_dataTerm[order(mirs_dataTerm$Freq, decreasing = TRUE),]
mirs_dataTerm <- subset(mirs_dataTerm, mirs_dataTerm[ , 2] > 2)  

# Get the miRNAs associated with the top GO terms
mir <- as.data.frame(mir)
mir <- mir[rep(seq_len(nrow(mir)), each = 10), ]
mirs_dataTerm$mirs <- 1
for (i in 1:nrow(mirs_dataTerm)) {
  term <- toString(mirs_dataTerm[i,1])
  ind <- which(mirs_data$Term == term)
  mirs_dataTerm[i,3] <- toString(mir[ind])
}
rm(i, ind, term)
write.table(mirs_dataTerm, "mirs_dataTerm.txt", append = FALSE, sep = "\t", dec = ".", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

# Get genes associated with pathways for each miRNA
mir <- read_excel("miRNAs_test.xlsx", col_names = FALSE)
colnames(mir) <- c("miRNAs")
mir <- mir$miRNAs
mirs_data2 <- mirs_data[c(1,2)]
mirs_data_merge <- merge(mirs_dataTerm, mirs_data2, by = "Term", sort = FALSE)
mirs_data_merge <- unique(mirs_data_merge)
go <- as.vector(mirs_data_merge$GO.ID)

k <- 1
for (i in 1:length(mir)){
  print(k)
  print(mir[i])
  predictions <- getPredictedTargets(mir[i], species = 'hsa', method = 'geom', min_src = 2)
  rankedGenes <- predictions[,'rank_product']
  selection <- function(x) TRUE
  allGO2genes <- annFUN.org(whichOnto='BP', feasibleGenes = NULL, mapping="org.Hs.eg.db", ID = "entrez")
  GOdata <- new('topGOdata', ontology = 'BP', allGenes = rankedGenes,
                annot = annFUN.GO2genes, GO2genes = allGO2genes, geneSel = selection, nodeSize=5)
  ann.genes <- genesInTerm(GOdata, go)
  output <- list()
  for(l in ann.genes){
    x <- getSYMBOL(l, data='org.Hs.eg')
    output <- append(output, list(x))
  }
  output <- as.data.frame(sapply(output, paste0, collapse=" "))
  rownames(output) <- NULL
  if (vector.is.empty(ann.genes)) {
    ann.genes <- c("0")
  }
  output$GO <- names(ann.genes)
  assign(paste(k, mir[i], sep = "_"), output)
  k <- k + 1
}
k <- 1

v <- paste(mirs_dataTerm$mirs, collapse=",")
v <- unique(trimws(unlist(strsplit(v, ","))))
for (i in 1:length(mir)){
  if (!mir[i] %in% v) {
    rm(list=ls(pattern = paste(mir[i], "$", sep = "")))
  }
}
for (i in ls(pattern = paste("^[1-9]"))){
  print(i)
  x <- get(i)
  colnames(x) <- c(paste(i), "GO.ID")
  assign(i, x)
}

list <- (ls(pattern = paste("^[1-9]")))
list[[length(list) + 1]] <- "mirs_data_merge"
multi_full <- as.data.frame(mirs_data_merge$GO.ID)
colnames(multi_full) <- c("GO.ID")

for (i in list){
  multi_full <- merge(multi_full, get(i), by = "GO.ID", all = TRUE)
}
for (i in 1:length(multi_full)){
  if (!mir[i] %in% v) {
    rm(list=ls(pattern = paste(mir[i], "$", sep = "")))
  }
}

colnames(multi_full)<-sub("^[0-9]*_","",colnames(multi_full))
# Select only miRNAs
df <- multi_full[2:21]

for (i in 1:nrow(df)){
  x <- colnames(df[,(names(df) %in% trimws(unlist(strsplit(paste(multi_full[i,ncol(multi_full)], collapse=","), ","))))])
  df[i, (which(!(colnames(df) %in% x)))] <- NA
}
df$GO.ID <- multi_full$GO.ID
df$Term <- multi_full$Term
df$mirs <- multi_full$mirs
write.table(df, "mirs_data_genes.txt", append = FALSE, sep = "\t", dec = ".", quote = FALSE,
                        row.names = FALSE, col.names = TRUE)
