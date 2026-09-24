#Load phenotypic data and SpaTS package####
library(SpATS)
library(dplyr)
library(data.table)
library(qqman)
library(ggplot2)
library(genio)


phen <- read.table("tlsCurveRatioBLUP.txt",sep="\t",header=T,as.is=T)
phen$VOL_3M <- phen$VOL_3M*0.01


#Create an empty matrix with rown the number of samples and columns the number of loops + 1
PhenoList <- matrix(NA, nrow = 774, ncol = 101)
#Names the first column <Trait> for Tassel and Pheno1, pheno2, pheno3 ... etc for the rest
colnames(PhenoList) <- c("<Trait>", paste("Pheno", 1:100, sep = ""))

#loops from 1 to the number of repeats
for (i in 1:100) {
  
  #So that the seed number changes each loop, otherwise all columns will be identical
  num = 68 * i
  #Seed is necesessery in order for others to be able to repeat the experiment.
  set.seed(num)
  
  #I make a new dataset from my pheno file that group all the clones together and then 
  #randomly select one of the clones each loop and put into the new dataset
  rando <- phen %>%
    group_by(Trait) %>%
    sample_n(1, replace = FALSE)
  
  #Put the random trees phenotype data into the phenoi column
  PhenoList[, i+1] <- rando$BLUPtlsCurveRatio
}

#In the first column <Trait>, list the Trait IDs
PhenoList[, 1] <- rando$Trait
#Now export this into a .txt file
write.table(PhenoList, "phenoBLUPCurveRatio.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#IN TERMINAL, DONT FORGET TO CHANGE PATH TO YOUR FILES

/Applications/TASSEL5/run_pipeline.pl -Xmx10g -log /Users/aanm0006/Documents/GWASLoop/BLUPloop/bootlog -fork1 -vcf /Users/aanm0006/Documents/GWASLoop/savar_old_only_biallelic_maf01_missing10_ordered.vcf -fork2 -r /Users/aanm0006/Documents/GWASLoop/phenoBLUPCurveRatio.txt -fork3 -k /Users/aanm0006/Documents/GWASLoop/kinshipmatrix.txt -combine4 -input1 -input2 -intersect -combine5 -input4 -input3 -mlm -mlmVarCompEst P3D -mlmCompressionLevel None -mlmOutputFile /Users/aanm0006/Documents/GWASLoop/BLUPloop/BLUPCurveRatioGWAS




