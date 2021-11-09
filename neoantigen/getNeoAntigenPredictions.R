library(biomaRt)
library(stringr)
library(ggplot2)
library(cowplot)

listEnsembl()
listEnsemblArchives()
listEnsembl(version = 102)
#listMarts()

fastaOutputWT="WT_proteinSequences9mers.fasta"
fastaOutputMUT="MUT_proteinSequences9mers.fasta"
mutMapName="MUT_proteinSequences9mers.mapFile"
wtMapName="WT_proteinSequences9mers.mapFile"

inputDir="/path/to/input/dir/"
variantAnnotations = "WES_RNA_merged_variant_annotations_filtered_missense_only.txt"

fastaFileWT<-paste(inputDir, fastaOutputWT, sep="")
fastaFileMUT<-paste(inputDir, fastaOutputMUT, sep="")
mutMap<-paste(inputDir, mutMapName, sep="")
wtMap<-paste(inputDir, wtMapName, sep="")

muts <- read.delim(paste(inputDir, variantAnnotations, sep = ""))
muts$Gene <- as.character(muts$Gene)
muts$Feature <- as.character(muts$Feature)
muts$Amino_acids <- as.character(muts$Amino_acids)
muts$X.Uploaded_variation <- as.character(muts$X.Uploaded_variation)

# MHC class I bind 9-11 aa peptides (run in WT and Mut and submit both to NetMHC so can compare differences in binding from WT and MUT)
kmerLength=9
window=kmerLength-1

hashListWT <- list()
hashListMUT <- list()
hashMapWT <- list()
hashMapMUT <- list()

ensembl=useMart("ensembl", host ="http://nov2020.archive.ensembl.org")
ensembl = useDataset("mmusculus_gene_ensembl", mart = ensembl)
attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)
listFilters(ensembl)
#getBM(attributes = c("ensembl_transcript_id","ensembl_gene_id","external_gene_name","mgi_symbol"),  filters = "ensembl_transcript_id", values = "ENSMUST00000044369", mart = ensembl)


# for (row in 1:nrow(muts)) {
#   protein <- getSequence(id = muts[row,"Feature"], type="ensembl_transcript_id", seqType = "peptide", mart=ensembl, verbose = FALSE)
#   kmerTargetRegion<-substring(protein$peptide, muts[row,"Protein_position"]-window, muts[row,"Protein_position"]+window) #Protein_position
#   if (paste(">",muts[row,"X.Uploaded_variation"],sep="") %in% names(hashListWT)) {
#     if (hashListWT[[ paste(">",muts[row,"X.Uploaded_variation"],sep="") ]]== kmerTargetRegion) {
#       print("already seen WT")
#     }else {
#       hashListWT[[ paste(">p",muts[row,"Protein_position"],"_",muts[row,"X.Uploaded_variation"], "_", muts[row,"Feature"], sep="") ]] <- kmerTargetRegion
#     }
#   }else {
#       hashListWT[[ paste(">",muts[row,"X.Uploaded_variation"],sep="") ]] <- kmerTargetRegion
#   }
#   
#   if (paste(">",muts[row,"X.Uploaded_variation"],sep="") %in% names(hashListMUT)) {
#     substr(kmerTargetRegion, kmerLength, kmerLength) <- strsplit(muts[row,"Amino_acids"], split = "/")[[1]][2]
#     if (hashListMUT[[ paste(">",muts[row,"X.Uploaded_variation"],sep="") ]]== kmerTargetRegion) {
#       print("already seen MUT")
#     }else {
#       hashListMUT[[ paste(">p",muts[row,"Protein_position"],"_",muts[row,"X.Uploaded_variation"], "_",muts[row,"Feature"], sep="") ]] <- kmerTargetRegion
#     }
#   }else {
#     substr(kmerTargetRegion, kmerLength, kmerLength) <- strsplit(muts[row,"Amino_acids"], split = "/")[[1]][2]
#     hashListMUT[[ paste(">",muts[row,"X.Uploaded_variation"],sep="") ]] <- kmerTargetRegion
#   }
# }

# testing start
#muts = muts[which(muts$Feature == "ENSMUST00000027810" | muts$Feature == "ENSMUST00000176740"),]

#wt1 <- getSequence(id = "ENSMUST00000027810", type="ensembl_transcript_id", seqType = "peptide", mart=ensembl, verbose = FALSE)
#wt2 <- getSequence(id = "ENSMUST00000176740", type="ensembl_transcript_id", seqType = "peptide", mart=ensembl, verbose = FALSE)

#kmer1<-substring(wt1, 296-8, 296+8) #Protein_position
#kmer2<-substring(wt2, 24-8, 24+8) #Protein_position


for (row in 1:nrow(muts)) {
  protein <- getSequence(id = muts[row,"Feature"], type="ensembl_transcript_id", seqType = "peptide", mart=ensembl, verbose = FALSE)
  kmerTargetRegion<-substring(protein$peptide, muts[row,"Protein_position"]-window, muts[row,"Protein_position"]+window) #Protein_position
  if (paste(muts[row,"X.Uploaded_variation"],"_",kmerTargetRegion,sep="") %in% names(hashMapWT)) {
      print("already seen WT")
      hashMapWT[[ paste(muts[row,"X.Uploaded_variation"],"_",kmerTargetRegion,sep="") ]] <- c(hashMapWT[[ paste(muts[row,"X.Uploaded_variation"],"_",kmerTargetRegion,sep="") ]], muts[row,"Feature"])
  }else {
    hashMapWT[[ paste(muts[row,"X.Uploaded_variation"],"_",kmerTargetRegion,sep="") ]] <- c(row, muts[row,"Feature"])
    hashListWT[[ paste(">",row,sep="") ]] <- kmerTargetRegion
  }

  substr(kmerTargetRegion, kmerLength, kmerLength) <- strsplit(muts[row,"Amino_acids"], split = "/")[[1]][2]
  if (paste(muts[row,"X.Uploaded_variation"],"_",kmerTargetRegion,sep="") %in% names(hashMapMUT)) {
      print("already seen MUT")
      hashMapMUT[[ paste(muts[row,"X.Uploaded_variation"],"_",kmerTargetRegion,sep="") ]] <- c(hashMapMUT[[ paste(muts[row,"X.Uploaded_variation"],"_",kmerTargetRegion,sep="") ]], muts[row,"Feature"])
  }else {
      hashMapMUT[[ paste(muts[row,"X.Uploaded_variation"],"_",kmerTargetRegion,sep="") ]] <- c(row, muts[row,"Feature"])
      hashListMUT[[ paste(">",row,sep="") ]] <- kmerTargetRegion
  }
}

# testing end



# write hashList as fasta file - WT
for (wtVariants in names(hashListWT)){
  cat(paste(wtVariants, hashListWT[[wtVariants]], sep="\n"), file=fastaFileWT, sep = "\n", append = T)
}

# write map file
for (wtVariants in names(hashMapWT)){
  cat(paste(wtVariants, paste( unlist(hashMapWT[[wtVariants]]), collapse=','), sep=","), file=wtMap, sep = "\n", append = T)
}

# write hashList as fasta file - MUT
for (mutVariants in names(hashListMUT)){
  cat(paste(mutVariants, hashListMUT[[mutVariants]], sep="\n"), file=fastaFileMUT, sep = "\n", append = T)
}

# write map file
for (mutVariants in names(hashMapMUT)){
  cat(paste(mutVariants, paste( unlist(hashMapMUT[[mutVariants]]), collapse=','), sep=","), file=mutMap, sep = "'\n", append = T)
}


# MHC class II bind 13-25 aa peptides (run on WT and Mut?)
kmerLength=9
window=kmerLength-1



############# POST NEOANTIGEN PREDICTION ####################
#strong binder = rank<= 0.5%
#weak binder = rank>0.5% and rank<= 2%
#no binding = rank>2%


mut_h2db<-read.csv("projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/proteinPrediction/MUT_proteinSequences9mers_MHC_class_I_H2Db_allele.tsv", sep = "\t", skip = 1)
wt_h2db <- read.csv("projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/proteinPrediction/WT_proteinSequences9mers_MHC_class_I_H2Db_allele.tsv", sep="\t", skip = 1)
mut_h2kb<-read.csv("projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/proteinPrediction/MUT_proteinSequences9mers_MHC_class_I_H2Kb_allele.tsv", sep = "\t", skip = 1)
wt_h2kb <- read.csv("projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/proteinPrediction/WT_proteinSequences9mers_MHC_class_I_H2Kb_allele.tsv", sep="\t", skip = 1)

mut_h2db$mhc_allele <- "H-2-Db"
wt_h2db$mhc_allele <- "H-2-Db"
mut_h2kb$mhc_allele <- "H-2-Kb"
wt_h2kb$mhc_allele <- "H-2-Kb"

mut_h2db$bind_strength <- "none"
wt_h2db$bind_strength <- "none"
mut_h2kb$bind_strength <- "none"
wt_h2kb$bind_strength <- "none"

mut_h2db$bind_strength[which(mut_h2db$Rank <= 0.5)] <- "strong"
wt_h2db$bind_strength[which(wt_h2db$Rank <= 0.5)] <- "strong"
mut_h2kb$bind_strength[which(mut_h2kb$Rank <= 0.5)] <- "strong"
wt_h2kb$bind_strength[which(wt_h2kb$Rank <= 0.5)] <- "strong"


mut_h2db$bind_strength[which(mut_h2db$Rank > 0.5 & mut_h2db$Rank <= 2.0)] <- "weak"
wt_h2db$bind_strength[which(wt_h2db$Rank > 0.5 & wt_h2db$Rank <= 2.0)] <- "weak"
mut_h2kb$bind_strength[which(mut_h2kb$Rank > 0.5 & mut_h2kb$Rank <= 2.0)] <- "weak"
wt_h2kb$bind_strength[which(wt_h2kb$Rank > 0.5 & wt_h2kb$Rank <= 2.0)] <- "weak"


names(mut_h2db) <- c("Pos", "Peptide_mut", "proteinID", "nM_mut", "Rank_mut", "Core_mut", "H_Avg_Ranks_mut", "N_binders_mut", "mhc_allele_mut", "bind_strength_mut")
names(mut_h2kb) <- c("Pos", "Peptide_mut", "proteinID", "nM_mut", "Rank_mut", "Core_mut", "H_Avg_Ranks_mut", "N_binders_mut", "mhc_allele_mut", "bind_strength_mut")
names(wt_h2db) <- c("Pos", "Peptide_wt", "proteinID", "nM_wt", "Rank_wt", "Core_wt", "H_Avg_Ranks_wt", "N_binders_wt", "mhc_allele_wt", "bind_strength_wt")
names(wt_h2kb) <- c("Pos", "Peptide_wt", "proteinID", "nM_wt", "Rank_wt", "Core_wt", "H_Avg_Ranks_wt", "N_binders_wt", "mhc_allele_wt", "bind_strength_wt")


mergePredsH2Db <- merge(x=mut_h2db, y = wt_h2db, by = c("proteinID", "Pos"), all = T)
mergePredsH2Kb <- merge(x=mut_h2kb, y = wt_h2kb, by = c("proteinID", "Pos"), all = T)

bindH2DbMutOnly <- mergePredsH2Db[which(mergePredsH2Db$N_binders_mut == 1 & mergePredsH2Db$N_binders_wt == 0),]
bindH2KbMutOnly <- mergePredsH2Kb[which(mergePredsH2Kb$N_binders_mut == 1 & mergePredsH2Kb$N_binders_wt == 0),]

bindH2DbBoth <- mergePredsH2Db[which(mergePredsH2Db$bind_strength_mut == "strong" & mergePredsH2Db$bind_strength_wt == "weak"),]
bindH2KbBoth <- mergePredsH2Kb[which(mergePredsH2Kb$bind_strength_mut == "strong" & mergePredsH2Kb$bind_strength_wt == "weak"),]

H2Db_peptides <- rbind(bindH2DbBoth, bindH2DbMutOnly)
H2Kb_peptides <- rbind(bindH2KbBoth, bindH2KbMutOnly)

map = read.csv("projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/proteinPrediction/map_dictionary.tsv", sep="\t")

H2Db_merged <- merge(H2Db_peptides, map, by = "proteinID", all.x = T)
H2Kb_merged <- merge(H2Kb_peptides, map, by = "proteinID", all.x = T)

# sample A223 txn counts -- move to new Test

############################## NEW TEST #########################################

# code for one sample of one allele

library(tidyverse)
txn <- read.csv("/home/tonya/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/quants/D4-RNA_quants/D4_quant.sf", sep="\t", header = T, col.names = c('Name', 'Length', 'EffectiveLength', 'D4_TPM', 'D4_NumReads'), colClasses = c("character"))
test <- as_tibble(txn) %>% separate(Name, into = c("txn", "subTxn"), "\\.", extra = "merge")
H2Db_merged$D4_counts <- ''
H2Db_merged$D4_txn_thresh <- 0
H2Db_merged$D4_txn_max <- 0
TPMthreshold = 1

for (row in 1:nrow(H2Db_merged)) {
  counts = ''
  getMaxTPM = c()
  search <- stringr::str_split(H2Db_merged[row, "Transcripts"], ";")[[1]]
  tmp <- test[which(test$txn %in% search),]
  for (i in seq(length(search))){
    getMaxTPM <- c(getMaxTPM, as.numeric(as.character(tmp$D4_TPM[i])))
    counts = paste(paste(tmp$txn[i],tmp$D4_TPM[i],tmp$D4_NumReads[i],sep="="), ";", counts, sep="")
  }
  H2Db_merged[row, "D4_counts"] <- counts
  if (max(getMaxTPM) > TPMthreshold){
    H2Db_merged[row, "D4_txn_thresh"] <- 1
    H2Db_merged[row, "D4_txn_max"] <- max(getMaxTPM)
  }
}

read_in_variants <- read.csv("projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/D4/D4_DNA_RNA_variant_info.tsv", sep="\t")
read_in_variants$SNP <- ''
read_in_variants$SNP <- paste(read_in_variants$ref, read_in_variants$alt, sep="/")
names(H2Db_merged)[19] <- "chromosome"
names(H2Db_merged)[20] <- "position"
D4_final <- merge(H2Db_merged, read_in_variants, by = c("chromosome", "position", "SNP"), all.x = T)
D4_final$RNA_totalReads_mean=rowMeans(D4_final[,c("D4_tophat2_RNA_totalReads", "D4_star_RNA_totalReads", "D4_hisat2_RNA_totalReads")], na.rm=TRUE)
D4_final_cleaned_H2Db <- D4_final[which((D4_final$D4_txn_thresh == 1) & (D4_final$D4_WES_DNA_totalReads >= 10) & (D4_final$RNA_totalReads_mean >=10)),]
write.table(D4_final_cleaned_H2Db, "~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/D4/D4_H2Db_final_neoantigen_candidates.txt", sep="\t", quote = F, row.names = F)

# other allele

txn <- read.csv("/home/tonya/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/quants/D4-RNA_quants/D4_quant.sf", sep="\t", header = T, col.names = c('Name', 'Length', 'EffectiveLength', 'D4_TPM', 'D4_NumReads'), colClasses = c("character"))
test <- as_tibble(txn) %>% separate(Name, into = c("txn", "subTxn"), "\\.", extra = "merge")
H2Kb_merged$D4_counts <- ''
H2Kb_merged$D4_txn_thresh <- 0
TPMthreshold = 1
H2Kb_merged$D4_txn_max <- 0

for (row in 1:nrow(H2Kb_merged)) {
  counts = ''
  getMaxTPM = c()
  search <- stringr::str_split(H2Kb_merged[row, "Transcripts"], ";")[[1]]
  tmp <- test[which(test$txn %in% search),]
  for (i in seq(length(search))){
    getMaxTPM <- c(getMaxTPM, as.numeric(as.character(tmp$D4_TPM[i])))
    counts = paste(paste(tmp$txn[i],tmp$D4_TPM[i],tmp$D4_NumReads[i],sep="="), ";", counts, sep="")
  }
  H2Kb_merged[row, "D4_counts"] <- counts
  if (max(getMaxTPM) > TPMthreshold){
    H2Kb_merged[row, "D4_txn_thresh"] <- 1
    H2Kb_merged[row, "D4_txn_max"] <- max(getMaxTPM)
  }
}

read_in_variants <- read.csv("projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/D4/D4_DNA_RNA_variant_info.tsv", sep="\t")
read_in_variants$SNP <- ''
read_in_variants$SNP <- paste(read_in_variants$ref, read_in_variants$alt, sep="/")
names(H2Kb_merged)[19] <- "chromosome"
names(H2Kb_merged)[20] <- "position"
D4_final <- merge(H2Kb_merged, read_in_variants, by = c("chromosome", "position", "SNP"), all.x = T)
D4_final$RNA_totalReads_mean=rowMeans(D4_final[,c("D4_tophat2_RNA_totalReads", "D4_star_RNA_totalReads", "D4_hisat2_RNA_totalReads")], na.rm=TRUE)
D4_final_cleaned_H2Kb <- D4_final[which((D4_final$D4_txn_thresh == 1) & (D4_final$D4_WES_DNA_totalReads >= 10) & (D4_final$RNA_totalReads_mean >=10)),]
write.table(D4_final_cleaned_H2Kb, "~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/D4/D4_H2Kb_final_neoantigen_candidates.txt", sep="\t", quote = F, row.names = F)



# all data merge
A223_H2Db <- read.csv("projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/A223/A223_H2Db_final_neoantigen_candidates.txt", sep="\t")
C12_H2Db <- read.csv("projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/C12/C12_H2Db_final_neoantigen_candidates.txt", sep="\t")
H10_H2Db <- read.csv("projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/H10/H10_H2Db_final_neoantigen_candidates.txt", sep="\t")
D4_H2Db <- read.csv("projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/D4/D4_H2Db_final_neoantigen_candidates.txt", sep="\t")

A223_H2Db$uniqueID <- paste(A223_H2Db$chromosome, A223_H2Db$position, A223_H2Db$SNP, A223_H2Db$Peptide_mut, sep="_")
C12_H2Db$uniqueID <- paste(C12_H2Db$chromosome, C12_H2Db$position, C12_H2Db$SNP, C12_H2Db$Peptide_mut, sep="_")
H10_H2Db$uniqueID <- paste(H10_H2Db$chromosome, H10_H2Db$position, H10_H2Db$SNP, H10_H2Db$Peptide_mut, sep="_")
D4_H2Db$uniqueID <- paste(D4_H2Db$chromosome, D4_H2Db$position, D4_H2Db$SNP, D4_H2Db$Peptide_mut, sep="_")

venn.diagram(x = list(A223_H2Db$uniqueID, C12_H2Db$uniqueID, D4_H2Db$uniqueID, H10_H2Db$uniqueID), 
             category.names = c("A223", "C12", "D4", "H10"), 
             filename = "H2Db_neoantigen_shared_across_clones.png", 
             imagetype = "png",
             height = 1000,
             width = 1000,
             col=c("#440154ff", '#21908dff', '#cc6677', '#FFA23A'),
             fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#cc6677',0.3), alpha('#FFA23A',0.3)),
             cex = 0.5,
             fontfamily = "sans",
             cat.cex = 0.7,
             cat.default.pos = "outer",
             cat.col = c("#440154ff", '#21908dff', '#cc6677', '#FFA23A'),
             cat.fontfamily = "sans", 
             main = "MHC-I H2Db")



# all data merge
A223_H2Kb <- read.csv("projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/A223/A223_H2Kb_final_neoantigen_candidates.txt", sep="\t")
C12_H2Kb <- read.csv("projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/C12/C12_H2Kb_final_neoantigen_candidates.txt", sep="\t")
H10_H2Kb <- read.csv("projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/H10/H10_H2Kb_final_neoantigen_candidates.txt", sep="\t")
D4_H2Kb <- read.csv("projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/D4/D4_H2Kb_final_neoantigen_candidates.txt", sep="\t")

A223_H2Kb$uniqueID <- paste(A223_H2Kb$chromosome, A223_H2Kb$position, A223_H2Kb$SNP, A223_H2Kb$Peptide_mut, sep="_")
C12_H2Kb$uniqueID <- paste(C12_H2Kb$chromosome, C12_H2Kb$position, C12_H2Kb$SNP, C12_H2Kb$Peptide_mut, sep="_")
H10_H2Kb$uniqueID <- paste(H10_H2Kb$chromosome, H10_H2Kb$position, H10_H2Kb$SNP, H10_H2Kb$Peptide_mut, sep="_")
D4_H2Kb$uniqueID <- paste(D4_H2Kb$chromosome, D4_H2Kb$position, D4_H2Kb$SNP, D4_H2Kb$Peptide_mut, sep="_")
library(VennDiagram)
venn.diagram(x = list(A223_H2Kb$uniqueID, C12_H2Kb$uniqueID, D4_H2Kb$uniqueID, H10_H2Kb$uniqueID), 
             category.names = c("A223", "C12", "D4", "H10"), 
             filename = "H2Kb_neoantigen_shared_across_clones.png", 
             imagetype = "png",
             height = 1000,
             width = 1000,
             col=c("#0072b2", '#f0e442', '#cc79a7', '#e69f00'),
             fill = c(alpha("#0072b2",0.3), alpha('#f0e442',0.3), alpha('#cc79a7',0.3), alpha('#e69f00',0.3)),
             cex = 0.5,
             fontfamily = "sans",
             cat.cex = 0.7,
             cat.default.pos = "outer",
             cat.col = c("#0072b2", '#f0e442', '#cc79a7', '#e69f00'),
             cat.fontfamily = "sans",
             main = "MHC-I H2Kb")

# bind strength prediction
A223_H2Db$strength <- paste(A223_H2Db$bind_strength_wt, A223_H2Db$bind_strength_mut, sep="->")
A223_H2Db$sample <- "A223"
C12_H2Db$strength <- paste(C12_H2Db$bind_strength_wt, C12_H2Db$bind_strength_mut, sep="->")
C12_H2Db$sample <- "C12"
H10_H2Db$strength <- paste(H10_H2Db$bind_strength_wt, H10_H2Db$bind_strength_mut, sep="->")
H10_H2Db$sample <- "H10"
D4_H2Db$strength<- paste(D4_H2Db$bind_strength_wt, D4_H2Db$bind_strength_mut, sep="->")
D4_H2Db$sample <- "D4"

allSamples <- rbind(A223_H2Db[,c('sample', 'strength')], C12_H2Db[,c('sample', 'strength')], H10_H2Db[,c('sample', 'strength')], D4_H2Db[,c('sample', 'strength')])

mutStengthSummary <- allSamples %>% group_by(sample, strength) %>% tally()
library(ggplot2)
library(viridis)

ggplot(mutStengthSummary, aes(fill=strength, y=n, x=sample)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) +
  ggtitle("H2Db Predicted Protein Bind Strength Based on Missense Mutation") +
  xlab("") +  labs(fill = "Predicted Bind Strength (WT -> MUT)") + ylab("Total number of predicted neoantigen") +
  theme(legend.position = "bottom", legend.box = "vertical", text = element_text(size=20), plot.title = element_text(hjust = 0.5))

# bind strength other allele

A223_H2Kb$strength <- paste(A223_H2Kb$bind_strength_wt, A223_H2Kb$bind_strength_mut, sep="->")
A223_H2Kb$sample <- "A223"
C12_H2Kb$strength <- paste(C12_H2Kb$bind_strength_wt, C12_H2Kb$bind_strength_mut, sep="->")
C12_H2Kb$sample <- "C12"
H10_H2Kb$strength <- paste(H10_H2Kb$bind_strength_wt, H10_H2Kb$bind_strength_mut, sep="->")
H10_H2Kb$sample <- "H10"
D4_H2Kb$strength<- paste(D4_H2Kb$bind_strength_wt, D4_H2Kb$bind_strength_mut, sep="->")
D4_H2Kb$sample <- "D4"

allSamples <- rbind(A223_H2Kb[,c('sample', 'strength')], C12_H2Kb[,c('sample', 'strength')], H10_H2Kb[,c('sample', 'strength')], D4_H2Kb[,c('sample', 'strength')])

mutStengthSummary <- allSamples %>% group_by(sample, strength) %>% tally()
library(ggplot2)
library(viridis)

ggplot(mutStengthSummary, aes(fill=strength, y=n, x=sample)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) +
  ggtitle("H2Kb Predicted Protein Bind Strength Based on Missense Mutation") +
  xlab("") +  labs(fill = "Predicted Bind Strength (WT -> MUT)") + ylab("Total number of predicted neoantigen") +
  theme(legend.position = "bottom", legend.box = "vertical", text = element_text(size=20), plot.title = element_text(hjust = 0.5))



# gene location of neoantigens
A223_H2Db$gene <- ''
for (row in 1:nrow(A223_H2Db)) {
  A223_H2Db[row, 'gene'] <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),  filters = "ensembl_transcript_id", values = str_split(A223_H2Db$Transcripts[row], ";")[[1]][1], mart = ensembl)$mgi_symbol
}

A223_H2Kb$gene <- ''
for (row in 1:nrow(A223_H2Kb)) {
  A223_H2Kb[row, 'gene'] <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),  filters = "ensembl_transcript_id", values = str_split(A223_H2Kb$Transcripts[row], ";")[[1]][1], mart = ensembl)$mgi_symbol
}

write.table(A223_H2Db, "~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/A223/A223_H2Db_final_neoantigen_candidates.txt", sep="\t", quote = F, row.names = F)
write.table(A223_H2Kb, "~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/A223/A223_H2Kb_final_neoantigen_candidates.txt", sep="\t", quote = F, row.names = F)



A223_H2Db_new <- A223_H2Db %>% group_by(strength, gene) %>% tally()
A223_H2Kb_new <- A223_H2Kb %>% group_by(strength, gene) %>% tally()

A223_H2Db_new$Sample <- ''
A223_H2Kb_new$Sample <- ''

A223_H2Db_new$Sample <- "A223"
A223_H2Kb_new$Sample <- "A223"

C12_H2Db$gene <- ''
for (row in 1:nrow(C12_H2Db)) {
  C12_H2Db[row, 'gene'] <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),  filters = "ensembl_transcript_id", values = str_split(C12_H2Db$Transcripts[row], ";")[[1]][1], mart = ensembl)$mgi_symbol
}

C12_H2Kb$gene <- ''
for (row in 1:nrow(C12_H2Kb)) {
  C12_H2Kb[row, 'gene'] <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),  filters = "ensembl_transcript_id", values = str_split(C12_H2Kb$Transcripts[row], ";")[[1]][1], mart = ensembl)$mgi_symbol
}
write.table(C12_H2Db, "~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/C12/C12_H2Db_final_neoantigen_candidates.txt", sep="\t", quote = F, row.names = F)
write.table(C12_H2Kb, "~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/C12/C12_H2Kb_final_neoantigen_candidates.txt", sep="\t", quote = F, row.names = F)

C12_H2Db_new <- C12_H2Db %>% group_by(strength, gene) %>% tally()
C12_H2Kb_new <- C12_H2Kb %>% group_by(strength, gene) %>% tally()

C12_H2Db_new$Sample <- ''
C12_H2Kb_new$Sample <- ''

C12_H2Db_new$Sample <- "C12"
C12_H2Kb_new$Sample <- "C12"


H10_H2Db$gene <- ''
for (row in 1:nrow(H10_H2Db)) {
  H10_H2Db[row, 'gene'] <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),  filters = "ensembl_transcript_id", values = str_split(H10_H2Db$Transcripts[row], ";")[[1]][1], mart = ensembl)$mgi_symbol
}

H10_H2Kb$gene <- ''
for (row in 1:nrow(H10_H2Kb)) {
  H10_H2Kb[row, 'gene'] <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),  filters = "ensembl_transcript_id", values = str_split(H10_H2Kb$Transcripts[row], ";")[[1]][1], mart = ensembl)$mgi_symbol
}

write.table(H10_H2Db, "~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/H10/H10_H2Db_final_neoantigen_candidates.txt", sep="\t", quote = F, row.names = F)
write.table(H10_H2Kb, "~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/H10/H10_H2Kb_final_neoantigen_candidates.txt", sep="\t", quote = F, row.names = F)

H10_H2Db_new <- H10_H2Db %>% group_by(strength, gene) %>% tally()
H10_H2Kb_new <- H10_H2Kb %>% group_by(strength, gene) %>% tally()


H10_H2Db_new$Sample <- ''
H10_H2Kb_new$Sample <- ''

H10_H2Db_new$Sample <- "H10"
H10_H2Kb_new$Sample <- "H10"

D4_H2Db$gene <- ''
for (row in 1:nrow(D4_H2Db)) {
  D4_H2Db[row, 'gene'] <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),  filters = "ensembl_transcript_id", values = str_split(D4_H2Db$Transcripts[row], ";")[[1]][1], mart = ensembl)$mgi_symbol
}

D4_H2Kb$gene <- ''
for (row in 1:nrow(D4_H2Kb)) {
  D4_H2Kb[row, 'gene'] <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),  filters = "ensembl_transcript_id", values = str_split(D4_H2Kb$Transcripts[row], ";")[[1]][1], mart = ensembl)$mgi_symbol
}

write.table(D4_H2Db, "~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/D4/D4_H2Db_final_neoantigen_candidates.txt", sep="\t", quote = F, row.names = F)
write.table(D4_H2Kb, "~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/D4/D4_H2Kb_final_neoantigen_candidates.txt", sep="\t", quote = F, row.names = F)


D4_H2Db_new <- D4_H2Db %>% group_by(strength, gene) %>% tally()
D4_H2Kb_new <- D4_H2Kb %>% group_by(strength, gene) %>% tally()

D4_H2Db_new$Sample <- ''
D4_H2Kb_new$Sample <- ''

D4_H2Db_new$Sample <- "D4"
D4_H2Kb_new$Sample <- "D4"

# upset plot
library(ComplexUpset)
# gene
test <- rbind(A223_H2Db_new, C12_H2Db_new, H10_H2Db_new, D4_H2Db_new)
none2strong <- test[which(test$strength == "none->strong"),c('gene', 'Sample')]

library(reshape2)

test <-melt(data = none2strong, measure.vars = c("gene")) %>% binarize()
test$Sample <- ''
test$Sample[which(test$Sample__A223 == 1)] <- "A223"
test$Sample[which(test$Sample__C12 == 1)] <- "C12"
test$Sample[which(test$Sample__D4 == 1)] <- "D4"
test$Sample[which(test$Sample__H10 == 1)] <- "H10"

names(test)[5:18] <- c("Apobec3", "Cenpe", "Dbf4", "Dguok", "Evc2", "Fanci", "Heatr5a", "Hjurp", "Kif13a", "Lipo3", "Mtus1", "Rasip1", "Sfxn4", "Xrcc4")

genes = colnames(test)[5:18]
upset(
  test,
  genes,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE,
      mapping=aes(fill=Sample)
    ) + scale_fill_manual(values=c(
      'A223'='#E41A1C', 'C12'='#377EB8',
      'H10'='#4DAF4A', 'D4'='#FF7F00'
    ))
  ),
  width_ratio=0.1
)

# NOTE this is something where interactive would be nice
#A223
ggplot(A223_H2Db, aes(x=A223_WES_DNA_vaf, y=RNA_avgVAF, color=strength, size=A223_txn_max)) + geom_point() + theme(legend.position="top", plot.title = element_text(hjust = 0.5)) + 
  labs(title="A223 H-2-Db VAF of Predicted Neoantigen DNA & RNA", x="WES VAF", y = "RNA VAF")

ggplot(A223_H2Kb, aes(x=A223_WES_DNA_vaf, y=RNA_avgVAF, color=strength, size=A223_txn_max)) + geom_point() + theme(legend.position="top", plot.title = element_text(hjust = 0.5)) + 
  labs(title="A223 H-2-Kb VAF of Predicted Neoantigen DNA & RNA", x="WES VAF", y = "RNA VAF")

#C12
ggplot(C12_H2Db, aes(x=C12_WES_DNA_vaf, y=RNA_avgVAF, color=strength, size=C12_txn_max)) + geom_point() + theme(legend.position="top", plot.title = element_text(hjust = 0.5)) + 
  labs(title="C12 H-2-Db VAF of Predicted Neoantigen DNA & RNA", x="WES VAF", y = "RNA VAF")

ggplot(C12_H2Kb, aes(x=C12_WES_DNA_vaf, y=RNA_avgVAF, color=strength, size=C12_txn_max)) + geom_point() + theme(legend.position="top", plot.title = element_text(hjust = 0.5)) + 
  labs(title="C12 H-2-Kb VAF of Predicted Neoantigen DNA & RNA", x="WES VAF", y = "RNA VAF")

#H10
ggplot(H10_H2Db, aes(x=H10_WES_DNA_vaf, y=RNA_avgVAF, color=strength, size=H10_txn_max)) + geom_point() + theme(legend.position="top", plot.title = element_text(hjust = 0.5)) + 
  labs(title="H10 H-2-Db VAF of Predicted Neoantigen DNA & RNA", x="WES VAF", y = "RNA VAF")

ggplot(H10_H2Kb, aes(x=H10_WES_DNA_vaf, y=RNA_avgVAF, color=strength, size=H10_txn_max)) + geom_point() + theme(legend.position="top", plot.title = element_text(hjust = 0.5)) + 
  labs(title="H10 H-2-Kb VAF of Predicted Neoantigen DNA & RNA", x="WES VAF", y = "RNA VAF")

#D4
ggplot(D4_H2Db, aes(x=D4_WES_DNA_vaf, y=RNA_avgVAF, color=strength, size=D4_txn_max)) + geom_point() + theme(legend.position="top", plot.title = element_text(hjust = 0.5)) + 
  labs(title="D4 H-2-Db VAF of Predicted Neoantigen DNA & RNA", x="WES VAF", y = "RNA VAF")

ggplot(D4_H2Kb, aes(x=D4_WES_DNA_vaf, y=RNA_avgVAF, color=strength, size=D4_txn_max)) + geom_point() + theme(legend.position="top", plot.title = element_text(hjust = 0.5)) + 
  labs(title="D4 H-2-Kb VAF of Predicted Neoantigen DNA & RNA", x="WES VAF", y = "RNA VAF")


############################## new anaysis ideas ##############################
A223_H2Kb <- read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/A223/A223_H2Kb_final_neoantigen_candidates.txt", sep = "\t")
C12_H2Kb <- read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/C12/C12_H2Kb_final_neoantigen_candidates.txt", sep = "\t")
H10_H2Kb <- read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/H10/H10_H2Kb_final_neoantigen_candidates.txt", sep = "\t")
D4_H2Kb <- read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/data_merge/D4/D4_H2Kb_final_neoantigen_candidates.txt", sep = "\t")



# H-2-Db
# get intersection of all and vaf
# inAllSamples <- Reduce(intersect, list(A223_H2Db$uniqueID, C12_H2Db$uniqueID, D4_H2Db$uniqueID, H10_H2Db$uniqueID))
A223_h2db_intersect<-A223_H2Db[which(A223_H2Db$strength == "none->strong"),c("A223_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]
C12_h2db_intersect<-C12_H2Db[which(C12_H2Db$strength == "none->strong"),c("C12_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]
H10_h2db_intersect <- H10_H2Db[which(H10_H2Db$strength == "none->strong"),c("H10_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]
D4_h2db_intersectV<- D4_H2Db[which(D4_H2Db$strength == "none->strong"),c("D4_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]
# 
# A223_h2db_intersect <- A223_H2Db[which(inAllSamples %in% A223_H2Db$uniqueID),c("A223_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]
# C12_h2db_intersect <- C12_H2Db[which(inAllSamples %in% C12_H2Db$uniqueID),c("C12_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]
# H10_h2db_intersect <- H10_H2Db[which(inAllSamples %in% H10_H2Db$uniqueID),c("H10_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]
# D4_h2db_intersect <- D4_H2Db[which(inAllSamples %in% D4_H2Db$uniqueID),c("D4_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]


names(A223_h2db_intersect)[1] <- 'WES_DNA_vaf'
names(C12_h2db_intersect)[1] <- 'WES_DNA_vaf'
names(H10_h2db_intersect)[1] <- 'WES_DNA_vaf'
names(D4_h2db_intersect)[1] <- 'WES_DNA_vaf'

library(gridExtra)
library(grid)
library(wesanderson)
grouped <- rbind(A223_h2db_intersect, C12_h2db_intersect, H10_h2db_intersect, D4_h2db_intersect)
grouped$uniqueID <- as.character(grouped$uniqueID)
strong <- grouped[which(grouped$strength == "none->strong"),]
myPlots <- list()
for (peptide in unique(strong$uniqueID)) {
  tmp <- strong[which(strong$uniqueID == peptide),]
  tmp1<-tmp[,c("RNA_avgVAF", "sample", "uniqueID", "gene")]
  tmp1$biotype <- "RNA"
  tmp2<-tmp[,c("WES_DNA_vaf", "sample", "uniqueID", "gene")]
  tmp2$biotype <- "DNA"
  names(tmp1)[1] <- "VAF"
  names(tmp2)[1] <- "VAF"
  tmp<-rbind(tmp1, tmp2)
  p1 <- ggplot(tmp, aes(fill=sample, y=VAF, x=biotype)) + geom_bar(position="dodge", stat="identity") + scale_fill_manual(values=wes_palette(n=4, name="Darjeeling2")) + 
          labs(title = paste(tmp$gene[1],":", " ",tmp$uniqueID[1], sep="")) + theme(plot.title = element_text(hjust = 0.5))
  myPlots[[peptide]] <- p1
}

do.call(grid.arrange, c(myPlots, nrow = 4))
# par(mfrow = c(5, 5))
# #listOfGraphs <- list()
# chr1<-grouped[which(startsWith(grouped$uniqueID, "1_")),]
# for (peptide in chr1$uniqueID) {
#   tmp <- chr1[which(chr1$uniqueID == peptide),]
#   print(ggplot(tmp, aes(fill=sample, y=RNA_avgVAF, x=as.factor(uniqueID))) + geom_bar(position="dodge", stat="identity"))
#   #listOfGraphs <- c(listOfGraphs, newtmp)
# }
# do.call(grid.arrange, c(listOfGraphs, nrow = 5))


# H-2-Kb
# get intersection of all and vaf
# inAllSamples <- Reduce(intersect, list(A223_H2Kb$uniqueID, C12_H2Kb$uniqueID, D4_H2Kb$uniqueID, H10_H2Kb$uniqueID))
A223_h2Kb_intersect<-A223_H2Kb[which(A223_H2Kb$strength == "none->strong"),c("A223_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]
C12_h2Kb_intersect<-C12_H2Kb[which(C12_H2Kb$strength == "none->strong"),c("C12_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]
H10_h2Kb_intersect <- H10_H2Kb[which(H10_H2Kb$strength == "none->strong"),c("H10_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]
D4_h2Kb_intersect<- D4_H2Kb[which(D4_H2Kb$strength == "none->strong"),c("D4_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]
# 
# A223_h2Kb_intersect <- A223_H2Kb[which(inAllSamples %in% A223_H2Kb$uniqueID),c("A223_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]
# C12_h2Kb_intersect <- C12_H2Kb[which(inAllSamples %in% C12_H2Kb$uniqueID),c("C12_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]
# H10_h2Kb_intersect <- H10_H2Kb[which(inAllSamples %in% H10_H2Kb$uniqueID),c("H10_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]
# D4_h2Kb_intersect <- D4_H2Kb[which(inAllSamples %in% D4_H2Kb$uniqueID),c("D4_WES_DNA_vaf", "RNA_avgVAF", "sample", "strength", "uniqueID", "gene")]


names(A223_h2Kb_intersect)[1] <- 'WES_DNA_vaf'
names(C12_h2Kb_intersect)[1] <- 'WES_DNA_vaf'
names(H10_h2Kb_intersect)[1] <- 'WES_DNA_vaf'
names(D4_h2Kb_intersect)[1] <- 'WES_DNA_vaf'

library(gridExtra)
library(grid)
library(wesanderson)
grouped <- rbind(A223_h2Kb_intersect, C12_h2Kb_intersect, H10_h2Kb_intersect, D4_h2Kb_intersect)
grouped$uniqueID <- as.character(grouped$uniqueID)
strong <- grouped[which(grouped$strength == "none->strong"),]
myPlots <- list()
for (peptide in unique(strong$uniqueID)) {
  tmp <- strong[which(strong$uniqueID == peptide),]
  tmp1<-tmp[,c("RNA_avgVAF", "sample", "uniqueID", "gene")]
  tmp1$biotype <- "RNA"
  tmp2<-tmp[,c("WES_DNA_vaf", "sample", "uniqueID", "gene")]
  tmp2$biotype <- "DNA"
  names(tmp1)[1] <- "VAF"
  names(tmp2)[1] <- "VAF"
  tmp<-rbind(tmp1, tmp2)
  p1 <- ggplot(tmp, aes(fill=sample, y=VAF, x=biotype)) + geom_bar(position="dodge", stat="identity") + scale_fill_manual(values=wes_palette(n=4, name="Royal1")) + 
    labs(title = paste(tmp$gene[1],":", " ",tmp$uniqueID[1], sep="")) + theme(plot.title = element_text(hjust = 0.5))
  myPlots[[peptide]] <- p1
}

do.call(grid.arrange, c(myPlots, nrow = 3))



########## separate analysis for RNA variant calling ##########
library(VennDiagram)
A223_star = read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/STEP8_independent_filtered_variant_calling/star/A223_star_vars.txt", header = F)
A223_tophat = read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/STEP8_independent_filtered_variant_calling/tophat2/A223_tophat2_vars.txt", header = F)
A223_hisat = read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/STEP8_independent_filtered_variant_calling/hisat2/A223_hisat_vars.txt", header = F)
venn.diagram(x = list(A223_star$V1, A223_tophat$V1, A223_hisat$V1), 
             category.names = c("STAR", "TopHat2", "HISAT2"), 
             filename = "A223_RNA_variant_calling_comparison.png", 
             imagetype = "png",
             height = 1000,
             width = 1000,
             col=c("#0072b2", '#f0e442', '#cc79a7'),
             fill = c(alpha("#0072b2",0.3), alpha('#f0e442',0.3), alpha('#cc79a7',0.3)),
             cex = 0.4,
             fontfamily = "sans",
             cat.cex = 0.5,
             cat.default.pos = "outer",
             cat.col = c("#0072b2", '#f0e442', '#cc79a7'),
             cat.fontfamily = "sans",
             main = "Sample: A223 variants called",
             main.cex = 0.5,
             print.mode = c("raw", "percent"))


C12_star = read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/STEP8_independent_filtered_variant_calling/star/C12_star_vars.txt", header = F)
C12_tophat = read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/STEP8_independent_filtered_variant_calling/tophat2/C12_tophat2_vars.txt", header = F)
C12_hisat = read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/STEP8_independent_filtered_variant_calling/hisat2/C12_hisat_vars.txt", header = F)
venn.diagram(x = list(C12_star$V1, C12_tophat$V1, C12_hisat$V1), 
             category.names = c("STAR", "TopHat2", "HISAT2"), 
             filename = "C12_RNA_variant_calling_comparison.png", 
             imagetype = "png",
             height = 1000,
             width = 1000,
             col=c("#0072b2", '#f0e442', '#cc79a7'),
             fill = c(alpha("#0072b2",0.3), alpha('#f0e442',0.3), alpha('#cc79a7',0.3)),
             cex = 0.4,
             fontfamily = "sans",
             cat.cex = 0.5,
             cat.default.pos = "outer",
             cat.col = c("#0072b2", '#f0e442', '#cc79a7'),
             cat.fontfamily = "sans",
             main = "Sample: C12 variants called",
             main.cex = 0.5,
             print.mode = c("raw", "percent"))


D4_star = read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/STEP8_independent_filtered_variant_calling/star/D4_star_vars.txt", header = F)
D4_tophat = read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/STEP8_independent_filtered_variant_calling/tophat2/D4_tophat2_vars.txt", header = F)
D4_hisat = read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/STEP8_independent_filtered_variant_calling/hisat2/D4_hisat_vars.txt", header = F)
venn.diagram(x = list(D4_star$V1, D4_tophat$V1, D4_hisat$V1), 
             category.names = c("STAR", "TopHat2", "HISAT2"), 
             filename = "D4_RNA_variant_calling_comparison.png", 
             imagetype = "png",
             height = 1000,
             width = 1000,
             col=c("#0072b2", '#f0e442', '#cc79a7'),
             fill = c(alpha("#0072b2",0.3), alpha('#f0e442',0.3), alpha('#cc79a7',0.3)),
             cex = 0.4,
             fontfamily = "sans",
             cat.cex = 0.5,
             cat.default.pos = "outer",
             cat.col = c("#0072b2", '#f0e442', '#cc79a7'),
             cat.fontfamily = "sans",
             main = "Sample: D4 variants called",
             main.cex = 0.5,
             print.mode = c("raw", "percent"))


H10_star = read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/STEP8_independent_filtered_variant_calling/star/H10_star_vars.txt", header = F)
H10_tophat = read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/STEP8_independent_filtered_variant_calling/tophat2/H10_tophat2_vars.txt", header = F)
H10_hisat = read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/STEP8_independent_filtered_variant_calling/hisat2/H10_hisat_vars.txt", header = F)
venn.diagram(x = list(H10_star$V1, H10_tophat$V1, H10_hisat$V1), 
             category.names = c("STAR", "TopHat2", "HISAT2"), 
             filename = "H10_RNA_variant_calling_comparison.png", 
             imagetype = "png",
             height = 1000,
             width = 1000,
             col=c("#0072b2", '#f0e442', '#cc79a7'),
             fill = c(alpha("#0072b2",0.3), alpha('#f0e442',0.3), alpha('#cc79a7',0.3)),
             cex = 0.4,
             fontfamily = "sans",
             cat.cex = 0.5,
             cat.default.pos = "outer",
             cat.col = c("#0072b2", '#f0e442', '#cc79a7'),
             cat.fontfamily = "sans",
             main = "Sample: H10 variants called",
             main.cex = 0.5,
             print.mode = c("raw", "percent"))

overall <- read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/RNA_alignment_metrics.csv")
ggplot(overall, aes(fill=Sample, y=Overall_Alignment, x=Aligner)) + geom_bar(position="dodge", stat="identity") + scale_fill_manual(values=wes_palette(n=4, name="Royal2")) +
  geom_text(position = position_dodge2(width = 1, preserve = "single"), aes(y=Overall_Alignment+1, label=Overall_Alignment, hjust=0.5), angle=0) +
  theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14))  

picard <- read.csv("~/projects/JingRachel_WES_and_RNAseq_variantCalling_11122020/RNA/picard_RNA_alignment_metrics.csv")

ggplot(picard, aes(y = ALIGNER, x = PCT_R1_TRANSCRIPT_STRAND_READS, label = PCT_R1_TRANSCRIPT_STRAND_READS, fill = SAMPLE, colour = SAMPLE)) +
  geom_segment(aes(x = 0, y = ALIGNER, xend =  PCT_R1_TRANSCRIPT_STRAND_READS, yend = ALIGNER), color = "grey50", size = 0.75) +
  geom_point(size = 3) +
  facet_wrap(~SAMPLE) + theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14))


ggplot(picard, aes(y = ALIGNER, x = PCT_R2_TRANSCRIPT_STRAND_READS, label = PCT_R2_TRANSCRIPT_STRAND_READS, fill = SAMPLE, colour = SAMPLE)) +
  geom_segment(aes(x = 0, y = ALIGNER, xend =  PCT_R2_TRANSCRIPT_STRAND_READS, yend = ALIGNER), color = "grey50", size = 0.75) +
  geom_point(size = 3) +
  facet_wrap(~SAMPLE)+ theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14))

ggplot(picard, aes(fill=Sample, y=Overall_Alignment, x=Aligner)) + geom_bar(position="dodge", stat="identity") + scale_fill_manual(values=wes_palette(n=4, name="Royal2")) +
  geom_text(position = position_dodge2(width = 1, preserve = "single"), aes(y=Overall_Alignment+1, label=Overall_Alignment, hjust=0.5), angle=0)

ggplot(picard, aes(x=SAMPLE, y=PCT_MRNA_BASES, fill=SAMPLE)) + geom_boxplot() + theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14))

ggplot(picard, aes(x=SAMPLE, y=PCT_RIBOSOMAL_BASES, fill=SAMPLE)) + geom_boxplot() + theme(axis.text.x = element_text(face="bold", size=14), axis.text.y = element_text(face="bold", size=14))
