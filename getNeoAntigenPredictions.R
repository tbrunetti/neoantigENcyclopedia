library(biomaRt)
library(stringr)

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



#strong binder = rank<= 0.5%
#weak binder = rank>0.5% and rank<= 2%
#no binding = rank>2%


mut_h2db<-read.csv("Downloads/software/netMHC-4.0/wt_mhc1_H2Db_9mer_predictions.tsv", sep = "\t")
wt_h2db <- read.csv("Downloads/software/netMHC-4.0/mut_mhc1_H2Db_9mer_predictions.tsv", sep="\t")
mut_h2kb<-read.csv("Downloads/software/netMHC-4.0/wt_mhc1_H2Kb_9mer_predictions.tsv", sep = "\t")
wt_h2kb <- read.csv("Downloads/software/netMHC-4.0/mut_mhc1_H2Kb_9mer_predictions.tsv", sep="\t")

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


names(mut_h2db) <- c("Pos", "Peptide_mut", "ID", "nM_mut", "Rank_mut", "Core_mut", "H_Avg_Ranks_mut", "N_binders_mut", "mhc_allele_mut", "bind_strength_mut")
names(mut_h2kb) <- c("Pos", "Peptide_mut", "ID", "nM_mut", "Rank_mut", "Core_mut", "H_Avg_Ranks_mut", "N_binders_mut", "mhc_allele_mut", "bind_strength_mut")
names(wt_h2db) <- c("Pos", "Peptide_wt", "ID", "nM_wt", "Rank_wt", "Core_wt", "H_Avg_Ranks_wt", "N_binders_wt", "mhc_allele_wt", "bind_strength_wt")
names(wt_h2kb) <- c("Pos", "Peptide_wt", "ID", "nM_wt", "Rank_wt", "Core_wt", "H_Avg_Ranks_wt", "N_binders_wt", "mhc_allele_wt", "bind_strength_wt")


mergePredsH2Db <- merge(x=mut_h2db, y = wt_h2db, by = c("ID", "Pos"), all = T)
mergePredsH2Kb <- merge(x=mut_h2kb, y = wt_h2kb, by = c("ID", "Pos"), all = T)

bindH2DbMutOnly <- mergePredsH2Db[which(mergePredsH2Db$N_binders_mut == 1 & mergePredsH2Db$N_binders_wt == 0),]
bindH2KbMutOnly <- mergePredsH2Kb[which(mergePredsH2Kb$N_binders_mut == 1 & mergePredsH2Kb$N_binders_wt == 0),]

bindH2DbBoth <- mergePredsH2Db[which(mergePredsH2Db$bind_strength_mut == "strong" & mergePredsH2Db$bind_strength_wt == "weak"),]
bindH2KbBoth <- mergePredsH2Kb[which(mergePredsH2Kb$bind_strength_mut == "strong" & mergePredsH2Kb$bind_strength_wt == "weak"),]

H2Db_peptides <- rbind(bindH2DbBoth, bindH2DbMutOnly)
H2Kb_peptides <- rbind(bindH2KbBoth, bindH2KbMutOnly)

combineAll <- rbind(H2Db_peptides, H2Kb_peptides)
combineAll$otherTx <- ""

combineAll$ID <- as.character(combineAll$ID)
for (vars in 1:nrow(combineAll)){
  print(vars)
  if (startsWith(combineAll[vars,"ID"], prefix = "p")){
    tmp <- stringr::str_split(combineAll[vars,"ID"], pattern = "_", n=2)
    print(tmp)
    combineAll[vars, "otherTx"] <- tmp[[1]][1]
    combineAll[vars, "ID"] <- tmp[[1]][2]
  }
}

combineAll %>% separate(ID, c("col1", "col2"), "_")
