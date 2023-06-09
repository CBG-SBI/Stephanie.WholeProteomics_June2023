library(PhosR)
library(dplyr)
setwd("/mnt/8TB/Projects/Projects/Melinda_Stephanie/Stephanie.WholeProteomics_June2023/")
pgroups <- read.table(file = "Data/combined_protein_proteome_fragpipe.tsv", header = T, sep = "\t",quote='',stringsAsFactors = FALSE,comment.char="")
pgroups.Intensity <- pgroups[,grepl("*LFQ.Intensity*", colnames(pgroups))]
summary(duplicated(pgroups$Protein.ID))
rownames(pgroups.Intensity) <- paste0(pgroups$Protein.ID,"_",pgroups$Gene)
pgroups.Intensity[pgroups.Intensity=="0"]<-NA
# delete rows with all NAs
pgroups.Intensity <- pgroups.Intensity[rowSums(is.na(pgroups.Intensity)) != ncol(pgroups.Intensity), ]
colnames(pgroups.Intensity) <- sub(".MaxLFQ.Intensity", "", colnames(pgroups.Intensity))
colnames(pgroups.Intensity) <- sub("X221012kw_SMaherProteome", "", colnames(pgroups.Intensity))

# Experimental design data frame 
library(readxl)
experimental_design <- as.data.frame(read_excel("Data/SM_Proteome_Annotation File.xlsx"))
colnames(experimental_design)[1] <- "LABEL"
experimental_design$LABEL <- gsub(" ",".",experimental_design$LABEL)
experimental_design$LABEL <- gsub(".MaxLFQ.Intensity","",experimental_design$LABEL)
experimental_design$LABEL <- gsub("221012kw_SMaherProteome","",experimental_design$LABEL)
intersect(experimental_design$LABEL,colnames(pgroups.Intensity))
rownames(experimental_design) <- experimental_design$LABEL
test <- data.frame(experimental_design$LABEL,colnames(pgroups.Intensity))
# reordering according to the experimental_design
pgroups.Intensity <- pgroups.Intensity[,experimental_design$LABEL]
# should be TRUE
all(colnames(pgroups.Intensity ) == row.names(experimental_design))

#Diagnostic plots

### Missing values pattern
keep <- pgroups.Intensity[apply(pgroups.Intensity, 1, function(x) any(is.na(x))), ] # select only proteins with NA values
missval <- ifelse(is.na(keep), 0, 1)
pheatmap::pheatmap(missval,  main = "Missing values pattern", 
                   color = c("black", "red"), show_rownames = F, annotation_col = experimental_design[,c(2,3,4,5,8,9)], fontsize_col = 4,
                   filename = "Results/Missing_values pattern.jpg", width=16, height=9)

# Visualization Number of proteins per sample
library(ggplot2)
binary <- ifelse(is.na(pgroups.Intensity), 0, 1)
temp <- as.data.frame(colSums(binary))
temp$experiment <- rownames(temp)
keep2 <- pgroups.Intensity[apply(pgroups.Intensity, MARGIN = 1, function(x) !all(is.na(x) == TRUE)), ]

png("Results/Number_of_proteins_per_sample.png", units="in", width=8, height=4.5, res=300)
ggplot(temp, aes(x = experiment, y = `colSums(binary)`)) + geom_bar(stat = "identity") + 
  geom_hline(yintercept = nrow(keep2)) + labs(title = "Proteins per sample", 
                                              x = "", y = "Number of proteins") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=3))
dev.off()

#Delete samples
selection <- grep("NLF2_S1_A3_1_11597", experimental_design$LABEL)
experimental_design <- experimental_design[-c(147),]
pgroups.Intensity$NLF2_S1_A3_1_11597 <- NULL
# should be TRUE
all(colnames(pgroups.Intensity ) == row.names(experimental_design))

#####################
### Preprocessing ###
#####################
temp <- pgroups.Intensity
temp <- temp[rowSums(is.na(temp)) != ncol(temp), ]
temp <- log2(temp)
temp <- PhosR::selectGrps(as.matrix(temp), experimental_design$Condition, 0.8, n=1)
boxplot(temp)
temp <- PhosR::scImpute(temp, 0.5, experimental_design$Condition)
temp <- PhosR::tImpute(temp)
proteins <- temp
proteins[is.na(proteins)] <- 0
saveRDS(proteins, file = "Data/Intensities.rds")
saveRDS(experimental_design, file = "Data/experimental_design.rds")
