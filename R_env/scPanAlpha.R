########
## Imports
########
# BiocManager::install("MAST")
# BiocManager::install("scater")
# BiocManager::install("zellkonverter")
# BiocManager::install("LoomExperiment")
# BiocManager::install("scran")
# BiocManager::install("NMF")
# BiocManager::install("GSEABase")
# devtools::install_github("cellgeni/sceasy")
library(scran)
library(scater)
library(Seurat)
library(reticulate)
library(MAST)
library(SingleCellExperiment)
library(lme4)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(data.table)
library(GGally)
library(zellkonverter)

# BiocManager::install("GO.db")

# options(mc.cores =20)

# setwd("R_env/")

############
## T1D #####
############

sca <- readH5AD("../h5ad/adata_t1d.final.h5ad")
#library(zellkonverter)
#library(scRNAseq)
# writeH5AD(sca,"adata_test_t1d.sca.h5ad")

sca.X <- assay(sca,"X")
sca.X <- log1p(expm1(sca.X) * colData(sca)$size_factors)

colData(sca)$n_genes <- scale(colData(sca)$n_genes)

table(colData(sca)$cell.type)


#Create data subset for alpha
sca_alpha <- subset(sca, with(colData(sca), cell.type=='alpha'))

#Filter out non-expressed genes in the subsets
print("Dimensions before subsetting:")
print(dim(sca_alpha))
print("")

sca_alpha_filt = sca_alpha[rowSums(assay(sca_alpha)) != 0, ]

print("Dimensions after subsetting:")
print(dim(sca_alpha_filt))

sca_alpha_filt2 <- filterLowExpressedGenes(sca_alpha,threshold=0.2)

print("Dimensions after subsetting:")
print(dim(sca_alpha_filt2))
table(colData(sca_alpha_filt2)$HPAP_id)
sca_alpha_filt2.1 <- subset(sca_alpha_filt2,!HPAP_id%in%c("HPAP-032","HPAP-034","HPAP-028",
                                                          "HPAP-023","HPAP-036",
                                                          "HPAP-039"))
table(colData(sca_alpha_filt2.1)$Disease)
table(colData(sca_alpha_filt2.1)$HPAP_id)
# adjust levels
colData(sca_alpha_filt2.1)[colData(sca_alpha_filt2.1)$HPAP_id=="HPAP-056","Tech"] = "10X.3"
colData(sca_alpha_filt2.1)$Tech <- droplevels(colData(sca_alpha_filt2.1)$Tech)
colData(sca_alpha_filt2.1)$HPAP_id <- droplevels(colData(sca_alpha_filt2.1)$HPAP_id)
#levels(colData(sca_alpha_filt2.1)$Tech) <- c("Chromium","Chromium","Fluidigm")
print(dim(sca_alpha_filt2.1))
table(colData(sca_alpha_filt2.1)$Disease,colData(sca_alpha_filt2.1)$Tech)


# sca_alpha_filt2.1 <- subset(sca_alpha_filt2,!HPAP_id%in%c("HPAP-032","HPAP-024","HPAP-028","HPAP-039","HPAP-020","HPAP-055"))



getDEs <- function(zlm,lrt="DiseaseT1D"){
  summaryCond_alpha <- summary(zlm, doLRT=lrt)
  
  #summaryCond_alpha <- summary(zlmCond_alpha, doLRT='DiseaseT1D')
  summaryDt_alpha <- summaryCond_alpha$datatable
  
  result_alpha <- merge(summaryDt_alpha[contrast==lrt & component=='H',.(primerid, `Pr(>Chisq)`)], #P-vals
                        summaryDt_alpha[contrast==lrt & component=='logFC', .(primerid, coef)],
                        by='primerid') #logFC coefficients
  
  #Correct for multiple testing (FDR correction) and filtering
  result_alpha[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  alpha_de = result_alpha[result_alpha$FDR<0.05,, drop=F]
  alpha_de = alpha_de[order(alpha_de$FDR),]
  alpha_de
  return(alpha_de)}


getAllDEs <- function(zlm,lrt="DiseaseT1D"){
  summaryCond_alpha <- summary(zlm, doLRT=lrt)
  summaryDt_alpha <- summaryCond_alpha$datatable
  result_alpha <- merge(summaryDt_alpha[contrast==lrt & component=='H',.(primerid, `Pr(>Chisq)`)], #P-vals
                        summaryDt_alpha[contrast==lrt & component=='logFC', .(primerid, coef)],
                        by='primerid') #logFC coefficients
  #Correct for multiple testing (FDR correction) and filtering
  result_alpha[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  return(result_alpha)
}



colData(sca_alpha_filt2.1)
# model with more balanced dataset
zlmCond_alpha2.1 <- zlm(formula = ~ Disease + n_genes + (1| HPAP_id) + Tech + Sex + BMI + Race + Age,
                        #zlmCond_alpha2.2 <- zlm(formula = ~ Disease + n_genes + (1| HPAP_id) + Tech,
                        #zlmCond_alpha2.1 <- zlm(formula = ~ Disease + n_genes + (1| HPAP_id),
                        sca=sca_alpha_filt2.1,
                        method='glmer', ebayes=FALSE,fitArgsD=list(nAGQ=0),parallel=TRUE)


alpha_allDe <- getAllDEs(zlmCond_alpha2.1)
alpha_allDe$symbols <- as.character(rowData(sca_alpha_filt2.1[alpha_allDe$primerid,])$gene_symbols)
alpha_allDe[(alpha_allDe$FDR<=0.05 & abs(coef)>=0.5),"primerid"]


write.table(alpha_allDe,"allDeT1D.tsv")
alpha_de2.1 <- getDEs(zlmCond_alpha2.1)
alpha_de2.1$symbols <- as.character(rowData(sca_alpha_filt2.1[alpha_de2.1$primerid,])$gene_symbols)
hist(abs(alpha_de2.1$coef))
sum(abs(alpha_de2.1$coef)>=.5,na.rm=T) # number of DEGs after FC filters
write.table(alpha_de2.1[abs(alpha_de2.1$coef)>=.5],"/home/bosi/analisi/panalphat2d/DEGs_T1D.tsv")
sum(alpha_de2.1[abs(alpha_de2.1$coef)>=.5]$coef>0) # up-reg
sum(alpha_de2.1[abs(alpha_de2.1$coef)>=.5]$coef<0) # down-reg
hist(alpha_de2.1$coef,xlim=c(-0.2,0.2))
alpha_de2.1$symbols <- as.character(rowData(sca_alpha_filt2.1[alpha_de2.1$primerid,])$gene_symbols)
alpha_de2.1
summary(zlmCond_alpha2.1)
write.table(alpha_de2.1,"alpha_T1DvsND_DEGs.tsv",sep="\t")


############
## T2D #####
############

setwd("scPanAlpha/h5ad/")
sca_t2d <- readH5AD("../h5ad/adata_t2d.final.h5ad")
sca_t2d.X <- assay(sca_t2d,"X")


freq_expressed <- 0.2
FCTHRESHOLD <- log2(1.5)

colData(sca_t2d)$size_factors

scaRaw_t2d <- FromMatrix(
  log1p(expm1(sca_t2d.X) * colData(sca_t2d)$size_factors),
  colData(sca_t2d),
  rowData(sca_t2d))


names(colData(scaRaw_t2d))[names(colData(scaRaw_t2d)) == "Disease"] <- "condition"
names(colData(scaRaw_t2d))[names(colData(scaRaw_t2d)) == "n_counts"] <- "libSize"
names(colData(scaRaw_t2d))[names(colData(scaRaw_t2d)) == "n_genes"] <- "nGeneOn"

# aheatmap(assay(scaRaw[1:1000,]), labRow='', annCol=as.data.frame(colData(scaRaw)[,c('condition', 'ourfilter')]), distfun='spearman')
set.seed(123)
plotPCA(scaRaw_t2d)




sca_t2d.X <- log1p(expm1(sca_t2d.X) * colData(sca_t2d)$size_factors)
colData(sca_t2d)$n_genes = scale(colData(sca_t2d)$n_genes)

#Create data subset for alpha
sca_alpha_t2d <- subset(sca_t2d, with(colData(sca_t2d), cell.type=='alpha'))

#Filter out non-expressed genes in the subsets
print("Dimensions before subsetting:")
print(dim(sca_alpha_t2d))
print("")

sca_alpha_t2d_filt = sca_alpha_t2d[rowSums(assay(sca_alpha_t2d)) != 0, ]

print("Dimensions after subsetting:")
print(dim(sca_alpha_t2d_filt))

sca_alpha_t2d_filt2 <- filterLowExpressedGenes(sca_alpha_t2d,threshold=0.2)
sca_alpha_t2d_filt2["ENSG00000171552"]
print("Dimensions after subsetting:")

print(dim(sca_alpha_t2d_filt2))
sum(colData(sca_alpha_t2d_filt2)$Disease=="ND")
table(colData(sca_alpha_t2d_filt2)$HPAP_id,colData(sca_alpha_t2d_filt2)$Disease)
table(colData(sca_alpha_t2d_filt2)$HPAP_id)
table(colData(sca_alpha_t2d_filt2)$HPAP_id,colData(sca_alpha_t2d_filt2)$Disease,colData(sca_alpha_t2d_filt2)$Tech)
table(colData(sca_alpha_t2d_filt2)$HPAP_id)
sca_alpha_t2d_filt2.1 <- sca_alpha_t2d_filt2
sca_alpha_t2d_filt2.1 <- subset(sca_alpha_t2d_filt2,!HPAP_id%in%c("HPAP-001","HPAP-006","HPAP-007","HPAP-014"))
sca_alpha_t2d_filt2.2 <- subset(sca_alpha_t2d_filt2,!HPAP_id%in%c("HPAP-001","HPAP-006"))
# sca_alpha_filt2.1 <- subset(sca_alpha_filt2,!HPAP_id%in%c("HPAP-032","HPAP-028","HPAP-039","HPAP-020","HPAP-055"))
# sca_alpha_t2d_filt2.1 <- sca_alpha_t2d_filt2
print(dim(sca_alpha_t2d_filt2.1))


table(colData(sca_alpha_t2d_filt2.1)$Disease,colData(sca_alpha_t2d_filt2.1)$HPAP_id)


zlmCond_alpha_t2d_2.1 <- zlm(formula = ~ Disease + n_genes + Tech + (1| HPAP_id) + Sex + BMI + Race + Age,
                             #zlmCond_alpha_t2d_2.1 <- zlm(formula = ~ Disease + n_genes + Tech + (1| HPAP_id),
                             sca=sca_alpha_t2d_filt2.1,
                             method='glmer', ebayes=FALSE,fitArgsD=list(nAGQ=0),parallel=TRUE)

alpha_de_t2d_2.1 <- getDEs(zlmCond_alpha_t2d_2.1,"DiseaseT2D")
alpha_de_t2d_2.1
sum(abs(alpha_de_t2d_2.1$coef>=.5))
sum(abs(alpha_de_t2d_2.1$coef<=-.5))
hist(abs(alpha_de_t2d_2.1$coef))



zlmCond_alpha_t2d_2.2 <- zlm(formula = ~ Disease + n_genes + Tech + (1| HPAP_id) + Sex + BMI + Race + Age,
                             #zlmCond_alpha_t2d_2.1 <- zlm(formula = ~ Disease + n_genes + Tech + (1| HPAP_id),
                             sca=sca_alpha_t2d_filt2.2,
                             method='glmer', ebayes=FALSE,fitArgsD=list(nAGQ=0),parallel=TRUE)

alpha_de_t2d_2.2 <- getDEs(zlmCond_alpha_t2d_2.2,"DiseaseT2D")
alpha_de_t2d_2.1
alpha_allDe_t2d_2.2 <- getAllDEs(zlmCond_alpha_t2d_2.2,"DiseaseT2D")

sum(abs(alpha_de_t2d_2.2$coef)>=.5)
hist(abs(alpha_de_t2d_2.2$coef))

alpha_allDe_t2d_2.2$symbols <- as.character(rowData(sca_alpha_t2d_filt2.2[alpha_allDe_t2d_2.2$primerid,])$gene_symbols)
write.table(alpha_allDe_t2d_2.2,"allDeT2D.tsv")


########################
# GSEA T1D
########################

kegg <- "c2.cp.kegg.v7.2.symbols.gmt"
reactome <- "c2.cp.reactome.v7.2.symbols.gmt"
go.bp <- "c5.go.bp.v7.2.symbols.gmt"
go.mf <- "c5.go.mf.v7.2.symbols.gmt"
go.cc <- "c5.go.cc.v7.2.symbols.gmt"
hallmark <- "h.all.v7.2.symbols.gmt"

set.seed(123)
# for reproducibility, I computed and saved boostrap data, which is now loaded
# the lines below describe how I computed the bootstrap

# boots <- bootVcov1(zlmCond_alpha2.1, R = 50)
# saveRDS(boots, 'bootstraps.rds')

boots <- readRDS("bootstraps.rds")

gmt2sigEnrichedModules <- function(gmt,
                                   sca=sca_alpha_filt2.1,
                                   zlm=zlmCond_alpha2.1,
                                   term="DiseaseT1D",
                                   min_gene_in_module=5,
                                   signThreshold=0.01){
  gmt <- getGmt(gmt)
  gene_ids <- geneIds(gmt)
  sets_indices <- limma::ids2indices(gene_ids, mcols(sca)$gene_symbols) # indicate the correct col name for gene symbols
  # Only keep modules with at least min_gene_in_module
  sets_indices <- sets_indices[sapply(sets_indices, length) >= min_gene_in_module]
  gsea <- gseaAfterBoot(zlm, boots, sets_indices, CoefficientHypothesis(term)) 
  z_stat_comb <- summary(gsea)
  sigModules <- z_stat_comb[combined_adj<signThreshold]
  return(sigModules)
}

gmt2enrichedModulesAll <- function(gmt,
                                   sca=sca_alpha_filt2.1,
                                   zlm=zlmCond_alpha2.1,
                                   term="DiseaseT1D",
                                   min_gene_in_module=5){
  gmt <- getGmt(gmt)
  gene_ids <- geneIds(gmt)
  sets_indices <- limma::ids2indices(gene_ids, mcols(sca)$gene_symbols) # indicate the correct col name for gene symbols
  # Only keep modules with at least min_gene_in_module
  sets_indices <- sets_indices[sapply(sets_indices, length) >= min_gene_in_module]
  gsea <- gseaAfterBoot(zlm, boots, sets_indices, CoefficientHypothesis(term)) 
  z_stat_comb <- summary(gsea)
  sigModules <- z_stat_comb
  return(sigModules)
}

gsea2plot <- function(gsea){
  gseaTable <- melt(gsea[1:20,.(set, combined_Z)], id.vars='set')
  gseaTable$set <- factor(gseaTable$set, levels = rev(gseaTable$set))
  ggplot(gseaTable, aes(y=set, x=variable, fill=value))+geom_raster() + scale_fill_distiller(palette="PiYG")
}

kegg.gsea <- gmt2sigEnrichedModules(kegg)
reactome.gsea <- gmt2sigEnrichedModules(reactome)
go.bp.gsea <- gmt2sigEnrichedModules(go.bp)
go.mf.gsea <- gmt2sigEnrichedModules(go.mf)
go.cc.gsea <- gmt2sigEnrichedModules(go.cc)
hallmark.gsea <- gmt2sigEnrichedModules(hallmark)

kegg.gsea.t1d.all <- gmt2enrichedModulesAll(kegg)
reactome.gsea.t1d.all <- gmt2enrichedModulesAll(reactome)
go.bp.gsea.t1d.all <- gmt2enrichedModulesAll(go.bp)
go.mf.gsea.t1d.all <- gmt2enrichedModulesAll(go.mf)
go.cc.gsea.t1d.all <- gmt2enrichedModulesAll(go.cc)
hallmark.gsea.t1d.all <- gmt2enrichedModulesAll(hallmark)


write.table(kegg.gsea,"/home/bosi/analisi/panalphat2d/kegg.gsea.t1d.tsv",sep = "\t")
write.table(reactome.gsea,"/home/bosi/analisi/panalphat2d/reactome.gsea.t1d.tsv",sep = "\t")
write.table(go.bp.gsea,"/home/bosi/analisi/panalphat2d/go.bp.gsea.t1d.tsv",sep = "\t")
write.table(go.mf.gsea,"/home/bosi/analisi/panalphat2d/go.mf.gsea.t1d.tsv",sep = "\t")
write.table(go.cc.gsea,"/home/bosi/analisi/panalphat2d/go.cc.gsea.t1d.tsv",sep = "\t")
write.table(hallmark.gsea,"/home/bosi/analisi/panalphat2d/hallmark.gsea.t1d.tsv",sep = "\t")


########################
# GSEA T2D
########################


# for reproducibility, I computed and saved boostrap data, which is now loaded
# the lines below describe how I computed the bootstrap

# boots_t2d <- bootVcov1(zlmCond_alpha_t2d_2.2, R = 50)
# saveRDS(boots_t2d, 'bootstraps_t2d.rds')

boots_t2d <- readRDS("bootstraps_t2d.rds")

kegg <- "c2.cp.kegg.v7.2.symbols.gmt"
reactome <- "c2.cp.reactome.v7.2.symbols.gmt"
go.bp <- "c5.go.bp.v7.2.symbols.gmt"
go.mf <- "c5.go.mf.v7.2.symbols.gmt"
go.cc <- "c5.go.cc.v7.2.symbols.gmt"
hallmark <- "h.all.v7.2.symbols.gmt"

gmt2sigEnrichedModules <- function(gmt,
                                   boots=boots,
                                   sca=sca_alpha_t2d_filt2.2,
                                   zlm=zlmCond_alpha_t2d_2.2,
                                   term="DiseaseT2D",
                                   min_gene_in_module=5){
  gmt <- getGmt(gmt)
  gene_ids <- geneIds(gmt)
  sets_indices <- limma::ids2indices(gene_ids, mcols(sca)$gene_symbols) # indicate the correct col name for gene symbols
  # Only keep modules with at least min_gene_in_module
  sets_indices <- sets_indices[sapply(sets_indices, length) >= min_gene_in_module]
  gsea <- gseaAfterBoot(zlm, boots, sets_indices, CoefficientHypothesis(term)) 
  z_stat_comb <- summary(gsea)
  sigModules <- z_stat_comb[combined_adj<.01]
  return(sigModules)
}

gmt2enrichedModulesAll <- function(gmt,
                                   boots=boots_t2d,
                                   sca=sca_alpha_t2d_filt2.2,
                                   zlm=zlmCond_alpha_t2d_2.2,
                                   term="DiseaseT2D",
                                   min_gene_in_module=5){
  gmt <- getGmt(gmt)
  gene_ids <- geneIds(gmt)
  sets_indices <- limma::ids2indices(gene_ids, mcols(sca)$gene_symbols) # indicate the correct col name for gene symbols
  # Only keep modules with at least min_gene_in_module
  sets_indices <- sets_indices[sapply(sets_indices, length) >= min_gene_in_module]
  gsea <- gseaAfterBoot(zlm, boots, sets_indices, CoefficientHypothesis(term)) 
  z_stat_comb <- summary(gsea)
  sigModules <- z_stat_comb
  return(sigModules)
}

gmt2sigEnrichedModules_ <- function(gmt,
                                    sca=sca_alpha_t2d_filt2.2,
                                    zlm=zlmCond_alpha_t2d_2.2,
                                    term="DiseaseT2D",
                                    min_gene_in_module=5){gmt2sigEnrichedModules(gmt,boots_t2d,sca,zlm,term,5)}  


kegg.gsea.t2d <- gmt2sigEnrichedModules_(kegg)
reactome.gsea.t2d <- gmt2sigEnrichedModules_(reactome)
go.bp.gsea.t2d <- gmt2sigEnrichedModules_(go.bp)
go.mf.gsea.t2d <- gmt2sigEnrichedModules_(go.mf)
go.cc.gsea.t2d <- gmt2sigEnrichedModules_(go.cc)
hallmark.gsea.t2d <- gmt2sigEnrichedModules_(hallmark)

kegg.gsea.t2d.all <- gmt2enrichedModulesAll(kegg)
reactome.gsea.t2d.all <- gmt2enrichedModulesAll(reactome)
go.bp.gsea.t2d.all <- gmt2enrichedModulesAll(go.bp)
go.mf.gsea.t2d.all <- gmt2enrichedModulesAll(go.mf)
go.cc.gsea.t2d.all <- gmt2enrichedModulesAll(go.cc)
hallmark.gsea.t2d.all <- gmt2enrichedModulesAll(hallmark)

write.table(kegg.gsea.t2d,"kegg.gsea.t2d.tsv",sep = "\t")
write.table(reactome.gsea.t2d,"reactome.gsea.t2d.tsv",sep = "\t")
write.table(go.bp.gsea.t2d,"go.bp.gsea.t2d.tsv",sep = "\t")
write.table(go.mf.gsea.t2d,"go.mf.gsea.t2d.tsv",sep = "\t")
write.table(go.cc.gsea.t2d,"go.cc.gsea.t2d.tsv",sep = "\t")
write.table(hallmark.gsea.t2d,"hallmark.gsea.t2d.tsv",sep = "\t")

############################################    
# Plot of GSEA (KEGG and Reactome)
############################################   
library(ggplot2)
## T1D

t1d.gsea_all <- read.table("Supplementary Table 5 - T1D GSEA - All.tsv",sep="\t",header = T)

t1d.gsea_kegg <- t1d.gsea_all[t1d.gsea_all$dataset == "KEGG",]
#t1d.gsea_kegg <- t1d.gsea_kegg[rev(order(abs(t1d.gsea_kegg$disc_effect))),]
t1d.gsea_kegg_ <- rbind(head(t1d.gsea_kegg[t1d.gsea_kegg$Enrichment.Sign==1,],n=10),
                        head(t1d.gsea_kegg[t1d.gsea_kegg$Enrichment.Sign==-1,],n=10))
t1d.gsea_kegg_$set <- str_replace_all(str_replace(t1d.gsea_kegg_$set, "KEGG_", ""),"_"," ")
t1d.gsea_kegg_$set <- factor(t1d.gsea_kegg_$set,levels=rev(t1d.gsea_kegg_$set))

pdf("t1d.gsea.kegg.pdf")
ggplot(t1d.gsea_kegg_, aes(y=set, fill=combined_adj_log, x=combined_Z)) +
  geom_bar(color="black",stat="identity") + 
  scale_fill_gradientn(colours=rev(c("#4575b4","#91bfdb","#e0f3f8"))) +
  theme(text = element_text(size=6),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black")) + 
  theme(panel.background = element_rect(fill = "white"))+
  labs(fill = "-log10(adjusted P-value)") +
  ggtitle("Top KEGG enriched terms in T1D vs ND") +
  xlab("Enrichment score") + ylab("Pathways")
dev.off()

t1d.gsea_reactome <- t1d.gsea_all[t1d.gsea_all$dataset == "Reactome",]
t1d.gsea_reactome_ <- rbind(head(t1d.gsea_reactome[t1d.gsea_reactome$Enrichment.Sign==1,],n=10),
                            head(t1d.gsea_reactome[t1d.gsea_reactome$Enrichment.Sign==-1,],n=10))
t1d.gsea_reactome_$set <- str_replace_all(str_replace(t1d.gsea_reactome_$set, "REACTOME_", ""),"_"," ")
t1d.gsea_reactome_$set <- factor(t1d.gsea_reactome_$set,levels=rev(t1d.gsea_reactome_$set))

pdf("t1d.gsea.reactome.pdf")
ggplot(t1d.gsea_reactome_, aes(y=set, fill=combined_adj_log, x=combined_Z)) +
  geom_bar(color="black",stat="identity") + 
  scale_fill_gradientn(colours=rev(c("#4575b4","#91bfdb","#e0f3f8"))) +
  theme(text = element_text(size=6),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black")) + 
  theme(panel.background = element_rect(fill = "white"))+
  labs(fill = "-log10(adjusted P-value)") +
  ggtitle("Top Reactome enriched terms in T1D vs ND") +
  xlab("Enrichment score") + ylab("Pathways")
dev.off()

t1d.msig <- t1d.gsea_all[t1d.gsea_all$dataset == "mSIG DB Hallmark",]
t1d.msig_ <- rbind(head(t1d.msig[t1d.msig$Enrichment.Sign==1,],n=10),
                   head(t1d.msig[t1d.msig$Enrichment.Sign==-1,],n=10))
t1d.msig_$set <- str_replace_all(str_replace(t1d.msig_$set, "HALLMARK_", ""),"_"," ")
t1d.msig_$set <- factor(t1d.msig_$set,levels=rev(t1d.msig_$set))

pdf("t1d.gsea.msig.pdf")
ggplot(t1d.msig_, aes(y=set, fill=combined_adj_log, x=combined_Z)) +
  geom_bar(color="black",stat="identity") + 
  scale_fill_gradientn(colours=rev(c("#4575b4","#91bfdb","#e0f3f8"))) +
  theme(text = element_text(size=6),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black")) + 
  theme(panel.background = element_rect(fill = "white"))+
  labs(fill = "-log10(adjusted P-value)") +
  ggtitle("Top MSIGDB enriched terms in T1D vs ND") +
  xlab("Enrichment score") + ylab("Pathways")
dev.off()

t1d.gobp <- t1d.gsea_all[t1d.gsea_all$dataset == "GO-BP",]
t1d.gobp_ <- rbind(head(t1d.gobp[t1d.gobp$Enrichment.Sign==1,],n=10),
                   head(t1d.gobp[t1d.gobp$Enrichment.Sign==-1,],n=10))
t1d.gobp_$set <- str_replace_all(str_replace(t1d.gobp_$set, "GO_", ""),"_"," ")
t1d.gobp_$set <- factor(t1d.gobp_$set,levels=rev(t1d.gobp_$set))

pdf("t1d.gsea.gobp.pdf")
ggplot(t1d.gobp_, aes(y=set, fill=combined_adj_log, x=combined_Z)) +
  geom_bar(color="black",stat="identity") + 
  scale_fill_gradientn(colours=rev(c("#4575b4","#91bfdb","#e0f3f8"))) +
  theme(text = element_text(size=6),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black")) + 
  theme(panel.background = element_rect(fill = "white"))+
  labs(fill = "-log10(adjusted P-value)") +
  ggtitle("Top GO (BP) enriched terms in T1D vs ND") +
  xlab("Enrichment score") + ylab("Pathways")
dev.off()

# T2D
t2d.gsea_all <- read.table("Supplementary Table 7 - T2D GSEA - All.tsv",sep="\t",header = T)

t2d.gsea_kegg <- t2d.gsea_all[t2d.gsea_all$dataset == "KEGG",]
t2d.gsea_kegg_ <- rbind(head(t2d.gsea_kegg[t2d.gsea_kegg$Enrichment.Sign==1,],n=10),
                        head(t2d.gsea_kegg[t2d.gsea_kegg$Enrichment.Sign==-1,],n=10))
t2d.gsea_kegg_$set <- str_replace_all(str_replace(t2d.gsea_kegg_$set, "KEGG_", ""),"_"," ")
t2d.gsea_kegg_$set <- factor(t2d.gsea_kegg_$set,levels=rev(t2d.gsea_kegg_$set))


pdf("t2d.gsea.kegg.pdf")
ggplot(t2d.gsea_kegg_, aes(y=set, fill=combined_adj_log, x=combined_Z)) +
  geom_bar(color="black",stat="identity") + 
  scale_fill_gradientn(colours=c('#fee0d2','#fc9272','#de2d26')) +
  theme(text = element_text(size=6),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black")) + 
  theme(panel.background = element_rect(fill = "white"))+
  labs(fill = "-log10(adjusted P-value)") +
  ggtitle("Top KEGG enriched terms in T2D vs ND") +
  xlab("Enrichment score") + ylab("Pathways")
dev.off()

t2d.gsea_reactome <- t2d.gsea_all[t2d.gsea_all$dataset == "Reactome",]
t2d.gsea_reactome_ <- rbind(head(t2d.gsea_reactome[t2d.gsea_reactome$Enrichment.Sign==1,],n=10),
                            head(t2d.gsea_reactome[t2d.gsea_reactome$Enrichment.Sign==-1,],n=10))
t2d.gsea_reactome_$set <- str_replace_all(str_replace(t2d.gsea_reactome_$set, "REACTOME_", ""),"_"," ")
t2d.gsea_reactome_$set <- factor(t2d.gsea_reactome_$set,levels=rev(t2d.gsea_reactome_$set))

pdf("t2d.gsea.reactome.pdf")
ggplot(t2d.gsea_reactome_, aes(y=set, fill=combined_adj_log, x=combined_Z)) +
  geom_bar(color="black",stat="identity") + 
  scale_fill_gradientn(colours=c('#fee0d2','#fc9272','#de2d26')) +
  theme(text = element_text(size=6),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black")) + 
  theme(panel.background = element_rect(fill = "white"))+
  labs(fill = "-log10(adjusted P-value)") +
  ggtitle("Top Reactome enriched terms in T2D vs ND") +
  xlab("Enrichment score") + ylab("Pathways")
dev.off()

t2d.msig <- t2d.gsea_all[t2d.gsea_all$dataset == "mSIG DB Hallmark",]
t2d.msig_ <- rbind(head(t2d.msig[t2d.msig$Enrichment.Sign==1,],n=10),
                   head(t2d.msig[t2d.msig$Enrichment.Sign==-1,],n=10))
t2d.msig_$set <- str_replace_all(str_replace(t2d.msig_$set, "HALLMARK_", ""),"_"," ")
t2d.msig_$set <- factor(t2d.msig_$set,levels=rev(t2d.msig_$set))

pdf("t2d.gsea.msig.pdf")
ggplot(t2d.msig_, aes(y=set, fill=combined_adj_log, x=combined_Z)) +
  geom_bar(color="black",stat="identity") + 
  scale_fill_gradientn(colours=c('#fee0d2','#fc9272','#de2d26')) +
  theme(text = element_text(size=6),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black")) + 
  theme(panel.background = element_rect(fill = "white"))+
  labs(fill = "-log10(adjusted P-value)") +
  ggtitle("Top mSIG enriched terms in T2D vs ND") +
  xlab("Enrichment score") + ylab("Pathways")
dev.off()

t2d.gobp <- t2d.gsea_all[t2d.gsea_all$dataset == "GO-BP",]
t2d.gobp_ <- rbind(head(t2d.gobp[t2d.gobp$Enrichment.Sign==1,],n=10),
                   head(t2d.gobp[t2d.gobp$Enrichment.Sign==-1,],n=10))
t2d.gobp_$set <- str_replace_all(str_replace(t2d.gobp_$set, "GO_", ""),"_"," ")
t2d.gobp_$set <- factor(t2d.gobp_$set,levels=rev(t2d.gobp_$set))

pdf("t2d.gsea.gobp.pdf",width=1150,height=637)
ggplot(t2d.gobp_, aes(y=set, fill=combined_adj_log, x=combined_Z)) +
  geom_bar(color="black",stat="identity") + 
  scale_fill_gradientn(colours=c('#fee0d2','#fc9272','#de2d26')) +
  theme(text = element_text(size=12),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black")) + 
  theme(panel.background = element_rect(fill = "white"))+
  labs(fill = "-log10(adjusted P-value)") +
  ggtitle("Top GO (BP) enriched terms in T2D vs ND") +
  xlab("Enrichment score") + ylab("Pathways")
dev.off()

