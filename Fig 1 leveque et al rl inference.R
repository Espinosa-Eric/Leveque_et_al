library (dplyr)
library (psych)
library(data.table)
library(ggplot2)


## load curated database from Soumelis according to Hou et al 

## this database can be downloaded at https://github.com/soumelis-lab/ICELLNET

LR2=read.csv2("D:/RNAseq/analysis/ICELLNETdb_up.csv",header = T)


colnames(LR2)[1:5]= c("Ligand.1", "Ligand.2","Receptor.1","Receptor.2","Receptor.3")



# load DEG helped vs non-helped MC (FDR<0.01 and LFC>1)
DEG_h_nh=read.csv2("D:/RNAseq/analysis/signif_res_h_nh.RL.csv", header = T, row.names=1)


DEG_h_nh=DEG_h_nh[!is.na(DEG_h_nh$symb),]
# load FPKM of Help condition from RNAseq
FPKM=read.csv2("D:/RNAseq/essais/mLT_B_MC_FPKM.csv", header=T)

# thresholding : keep only  genes showing minimal expression set to FPKM>1 according to the distribution
p <- ggplot(FPKM, aes(x=log10(FPKM+0.001))) + 
  stat_density()
p

help_FPKM=FPKM[!FPKM$FPKM<1,]



# get FPKM for DEG
DEG_h_nh$FPKM<- help_FPKM$FPKM[match(DEG_h_nh$symbol, help_FPKM$genesymb)]


# feed the RL datatable with MC DEG as ligands and FPKM for activated CD4+ T cells
# Whenever ligand or receptor  is not expressed, the score of this particular interaction is set to zero
#  MC genes FPKM input as ligands for DEG present in the R-L dataset 
LR2$hMC.1<- DEG_h_nh$FPKM[match(LR2$Ligand.1,DEG_h_nh$symbol)]
LR2$hMC.1[is.na(LR2$hMC.1)] <- 0
LR2$hMC.2<- DEG_h_nh$FPKM[match(LR2$Ligand.2,DEG_h_nh$symbol)]
LR2$hMC.2[is.na(LR2$hMC.2)] <- 0

# indicate in the R-L datatable the log2FoldChange (calculated by Deseq2 for MC DEGs)
LR2$lfc<- DEG_h_nh$log2FoldChange[match(LR2$Ligand.1,DEG_h_nh$symbol)]
LR2$lfc[is.na(LR2$hMC.1)] <- 0




# Get the FPKM from CD4+ memory  T cell activated with  anti-CD3/28 beads for 24h
# from GSE73214
# PMID: 28947543

t1=read.table("D:/RNAseq/essais/GSM1888832_5291-Memory-T24.txt", header = T)
t2=read.table("D:/RNAseq/essais/GSM1888829_4659-Memory-T24.txt", header = T)
t3=read.table("D:/RNAseq/essais/GSM1888830_5053-Memory-T24.txt", header = T)
t4=read.table("D:/RNAseq/essais/GSM1888831_5131-Memory-T24.txt", header = T)


moyFPKM=(t1$FPKM+t2$FPKM+t3$FPKM+t4$FPKM)/4


act.CD4_FPKM=data.frame(cbind(t1[,2],moyFPKM))
act.CD4_FPKM$moyFPKM=as.numeric(act.CD4_FPKM$moyFPKM)
colnames(act.CD4_FPKM)[1]="Receptor.gene.symbol"



# thresholding : keep only  genes showing minimal expression set to FPKM>1

q <- ggplot(act.CD4_FPKM, aes(x=log10(moyFPKM+0.001))) + 
  stat_density()
q

act.CD4_FPKM=act.CD4_FPKM[!act.CD4_FPKM$moyFPKM<1,]


# Add in the R-L datatable the FPKM of matched Receptor in act.CD4_FPKM 

LR2$TCD4.1<- act.CD4_FPKM[match(LR2$Receptor.1,act.CD4_FPKM$Receptor.gene.symbol),2]
LR2$TCD4.1[is.na(LR2$TCD4.1)] <- 0

LR2$TCD4.2<- act.CD4_FPKM[match(LR2$Receptor.2,act.CD4_FPKM$Receptor.gene.symbol),2]
LR2$TCD4.2[is.na(LR2$TCD4.2)] <- 0

LR2$TCD4.3<- act.CD4_FPKM[match(LR2$Receptor.3,act.CD4_FPKM$Receptor.gene.symbol),2]
LR2$TCD4.3[is.na(LR2$TCD4.3)] <- 0


#After selecting the genes corresponding to the ligands and/or receptors from the transcriptional profiles, 
#each ligand or  receptor gene expression is scaled by a maximum of gene expression and then multiplied by 10, 
#to have values ranging from 0 to 10. For each gene, the maximum value (10) is defined as the mean of expression of the 10%
#highest values of expression for RNA-seq and microarray datasets. Outliers are
#rescaled at 10 if above maximum value.




FPKM_Lig=LR2$hMC.1[!LR2$hMC.1<=0]
upper_boundL <- quantile(FPKM_Lig, 0.9)
upper_boundL

outlier_indL <- which(LR2$hMC.1 > upper_boundL)
outlier_indL

maxL= mean(LR2$hMC.1[outlier_indL])

LR2$hMC.1=LR2$hMC.1/maxL*10
LR2$hMC.1[which(LR2$hMC.1>10)]=10

LR2$hMC.2=LR2$hMC.2/maxL*10
LR2$hMC.2[which(LR2$hMC.2>10)]=10

FPKM_Rec=LR2$TCD4.1[!LR2$TCD4.1<=0]
upper_boundR <- quantile(FPKM_Rec, 0.9)
upper_boundR

outlier_indR <- which(LR2$TCD4.1 > upper_boundR)
outlier_indR

maxR= mean(LR2$TCD4.1[outlier_indR])

LR2$TCD4.1=LR2$TCD4.1/maxR*10
LR2$TCD4.1[which(LR2$TCD4.1>10)]=10

LR2$TCD4.2=LR2$TCD4.2/maxR*10
LR2$TCD4.2[which(LR2$TCD4.2>10)]=10

LR2$TCD4.3=LR2$TCD4.3/maxR*10
LR2$TCD4.3[which(LR2$TCD4.3>10)]=10

# score calculation with geomean for multichain receptors


for (mol in seq(1, dim(LR2)[1])) {
     # In receptor-provider data
    if (is.na(LR2$Receptor.2[mol]) & is.na(LR2$Receptor.3[mol])) {
      LR2$Rscore[mol] = LR2$TCD4.1[mol]
    } else if (is.na(LR2$Receptor.3[mol])) {
      LR2$Rscore[mol] = geometric.mean(c(LR2$TCD4.1[mol], LR2$TCD4.2[mol]))
    } else if (is.na(LR2$Receptor.2[mol])) {
      LR2$Rscore[mol] = geometric.mean(c(LR2$TCD4.1[mol], LR2$TCD4.3[mol]))
    } else{
      LR2$Rscore[mol] = geometric.mean(c(LR2$TCD4.1[mol], LR2$TCD4.2[mol],LR2$TCD4.3[mol]) )}
}

# score for multichian ligands in helped MC

for (mol in seq(1, dim(LR2)[1])) {
  # In ligand-provider data
  if (is.na(LR2$Ligand.2[mol]) ) {
    LR2$hLscore[mol] = LR2$hMC.1[mol]
   } else{
    LR2$hLscore[mol] = geometric.mean(c(LR2$hMC.1[mol], LR2$hMC.2[mol]) )}
}


# global score 
 LR2$hRLscore=LR2$Rscore * LR2$hLscore


 
 # extract the R-L pairs having communication score>0
 
 res=LR2[which(LR2$hRLscore>0),]
 write.csv2(res,"D:/RNAseq/analysis/table_RL_scores.csv")
 # select the top 30 R-L score without HLA ligand that would overrepresent the others in the circosplot 
 res_to_plot=res[-grep("HLA",res$Ligand.1, value = F),]
 
 res_to_plot= res_to_plot %>% top_n(30,hRLscore)

 res_to_plot=res_to_plot[,c(1,3,25)]
 res_to_plot$hRLscore=scale.default(res_to_plot$hRLscore)

 library(DT)
 datatable(res)

# circosplot representation of the top 30 Receptor-Ligand pair inferred
library(circlize)


circos.par(start.degree = -90)
chordDiagram(res_to_plot, directional = 1,direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", diffHeight = mm_h(2),
             annotationTrack = "grid",
             
             preAllocateTracks = 1)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.7)
}, bg.border = NA)

circos.clear()


