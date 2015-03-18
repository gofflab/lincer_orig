#!/usr/bin/env Rscript --vanilla

########################
# This is by no means a working script or template.  
# Please adjust as needed for you arrays.
#
# Loyal A. Goff 2011 
########################

#Load libraries

library(limma)

#Load Targets
targets <- readTargets("Targets.txt")

#Remove blank sample
#targets<-targets[-16,]

#Get Dave's updated annotation
annot<-read.delim("Array_Annotation_tgs.txt",header=T,sep="\t",na.string="")

#Read Arrays
RG<-read.maimages(targets,columns=list(G="gMedianSignal",Gb = "gBGMedianSignal",R="rMedianSignal",Rb="rBGMedianSignal"),
					annotation=c("Row","Col","FeatureNum","ControlType","ProbeName","GeneName","SystematicName")
					)
					
#Apply Dave's Updated annotation #---Make sure annotation is in correct order or labels will be corrupted---#
RG$genes$GeneName<-annot$GeneName
RG$genes$SystematicName<-annot$SystematicName

#Background Subtraction
RG <- backgroundCorrect(RG, method="minimum", offset=1)

#NormalizeWithinArrays
MA <- normalizeWithinArrays(RG, method="loess")

#NormalizeBetweenArrays
MA <- normalizeBetweenArrays(MA, method="quantile")

#Average Replicates
MA.avg <- avereps(MA, ID=MA$genes$GeneName)

#########
#######
#LIMMA
######
#########


#####
#vs pTrex
#######

#Create Design Matrix
design <- modelMatrix(targets, ref="pTrex")

#lmFit
fit <- lmFit(MA.avg, design)
fit2 <- eBayes(fit)

#Test Tables
output <- topTable(fit2, adjust="BH", coef="CW1.1", number=40000)
write.table(output, file="CW1.1_vs_pTrex.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="CW1.1_vs_pTrex.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

output <- topTable(fit2, adjust="BH", coef="CW1.2", number=40000)
write.table(output, file="CW1.2_vs_pTrex.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="CW1.2_vs_pTrex.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

output <- topTable(fit2, adjust="BH", coef="CW1.3", number=40000)
write.table(output, file="CW1.3_vs_pTrex.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="CW1.3_vs_pTrex.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

output <- topTable(fit2, adjust="BH", coef="CW2.1", number=40000)
write.table(output, file="CW2.1_vs_pTrex.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="CW2.1_vs_pTrex.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

output <- topTable(fit2, adjust="BH", coef="CW3.1", number=40000)
write.table(output, file="CW3.1_vs_pTrex.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="CW3.1_vs_pTrex.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

output <- topTable(fit2, adjust="BH", coef="CW6.1", number=40000)
write.table(output, file="CW6.1_vs_pTrex.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="CW6.1_vs_pTrex.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

output <- topTable(fit2, adjust="BH", coef="UTC", number=40000)
write.table(output, file="UTC_vs_pTrex.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="UTC_vs_pTrex.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

output <- topTable(fit2, adjust="BH", coef="GFP", number=40000)
write.table(output, file="GFP_vs_pTrex.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="GFP_vs_pTrex.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

#####
#vs GFP
#######

#Create Design Matrix
design <- modelMatrix(targets, ref="GFP")

#lmFit
fit <- lmFit(MA.avg, design)
fit2 <- eBayes(fit)

#Test Tables
output <- topTable(fit2, adjust="BH", coef="CW1.1", number=40000)
write.table(output, file="CW1.1_vs_GFP.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="CW1.1_vs_GFP.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

output <- topTable(fit2, adjust="BH", coef="CW1.2", number=40000)
write.table(output, file="CW1.2_vs_GFP.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="CW1.2_vs_GFP.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

output <- topTable(fit2, adjust="BH", coef="CW1.3", number=40000)
write.table(output, file="CW1.3_vs_GFP.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="CW1.3_vs_GFP.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

output <- topTable(fit2, adjust="BH", coef="CW2.1", number=40000)
write.table(output, file="CW2.1_vs_GFP.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="CW2.1_vs_GFP.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

output <- topTable(fit2, adjust="BH", coef="CW3.1", number=40000)
write.table(output, file="CW3.1_vs_GFP.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="CW3.1_vs_GFP.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

output <- topTable(fit2, adjust="BH", coef="CW6.1", number=40000)
write.table(output, file="CW6.1_vs_GFP.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="CW6.1_vs_GFP.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

output <- topTable(fit2, adjust="BH", coef="UTC", number=40000)
write.table(output, file="UTC_vs_GFP.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="UTC_vs_GFP.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

output <- topTable(fit2, adjust="BH", coef="pTrex", number=40000)
write.table(output, file="pTrex_vs_GFP.txt", sep="\t", quote=FALSE)
write.table(output[,c(6,8)], file="pTrex_vs_GFP.rnk", row.names=F,col.names=F,sep="\t", quote=FALSE)

##############
#Output all with F-test for significance (vs T-test for individual coefs)
#output <- topTable(fit2, adjust="BH", coef=c("CW1.1","CW1.2","CW1.3","CW2.1","CW3.1","CW6.1","GFP"), number=40000, p.value=0.05)
#write.table(output, file="All_vs_pTrex.txt", sep="\t", quote=FALSE)
