# KT007_compare_array_to_ELISA

```r
rm(list=ls())
setwd("~/Desktop/Microarray Analyses/Kevin's array/151118 data/")
print(load("./Data/Processed Data/006.KTarray.Apac.clinical_and_array.RData"))
```

##############################################################################
##standardize naming structure for ELISA data
```r
names(apac_clinical)[grep("ama", names(apac_clinical))]
names(apac_clinical)[grep("msp1", names(apac_clinical))]
names(apac_clinical)[grep("msp2", names(apac_clinical))]
names(apac_clinical)[grep("csp", names(apac_clinical))]

names(apac_clinical) <- gsub("ama", "ama_", names(apac_clinical))
names(apac_clinical) <- gsub("ama_t_x1_n", "ama_t_x1", names(apac_clinical))
names(apac_clinical) <- gsub("msp1p", "msp1_p", names(apac_clinical))
names(apac_clinical) <- gsub("msp2_t_x", "msp2_pos_x", names(apac_clinical))
names(apac_clinical) <- gsub("msp2", "msp2_", names(apac_clinical))
names(apac_clinical) <- gsub("csp", "csp__", names(apac_clinical))
names(apac_clinical) <- gsub("[_][_]", "_", names(apac_clinical))
names(apac_clinical) <- gsub("__", "_", names(apac_clinical))
names(apac_clinical) <- gsub("apos", "ama_pos", names(apac_clinical))
names(apac_clinical) <- gsub("pos_t_", "pos_", names(apac_clinical))
names(apac_clinical) <- gsub("pos_t", "pos_x", names(apac_clinical))
```

##############################################################################
##find out array spots that correspond to ama1, msp1, msp2, and csp antigens
```r
ama_spots <- apac_spots[grep("AMA", toupper(apac_spots$description)),]
ama_spots <- ama_spots[-grep("GAMA", toupper(ama_spots$description)),]
msp1_spots <- apac_spots[grep("MSP1", toupper(apac_spots$description)),]
msp1_spots <- msp1_spots[-grep("MSP11", toupper(msp1_spots$description)),]
msp2_spots <- apac_spots[grep("MSP2", toupper(apac_spots$description)),]
csp_spots <- apac_spots[grep("CSP", toupper(apac_spots$description)),]

spots <- rbind(ama_spots, msp1_spots)
spots <- rbind(spots, msp2_spots)
spots <- rbind(spots, csp_spots)
```

##############################################################################
##subset out elisa titres, and set minimum to 0.1
```r
elisa <- apac_clinical[,c("sampleID", "study", "Sample", names(apac_clinical)[grep("_t_", names(apac_clinical))])]
for(i in colnames(elisa)[4:ncol(elisa)]){
  elisa[,i][elisa[,i]<=0 & !is.na(elisa[,i])] <- 0.1
}
```

##############################################################################
##merge eliza and array data
```r
expr <- intensity.MedA[rownames(intensity.MedA) %in% c(spots$Unique.ID),]
data <- data.frame(t(expr))
data$sampleID <- rownames(data)

elisa <- merge(elisa, data, by="sampleID")
```

##############################################################################
##make scatter plots
```r
pdf(file="./Figures/Exploratory/0007.KTarray.Apac_X1.scatter_elisa_vs_array.pdf",width=12,height=9)
par(mfcol=c(2,3))
x1_elisa <- elisa[elisa$study=="Apac_X1",]
for(antigen in c("ama", "msp1", "msp2", "csp")){
  microarray_spots <- spots[grep(toupper(antigen), spots$description),]
  for(ag in microarray_spots$UID){
    print(antigen)
    print(ag)
    plot(x=log(x1_elisa[,paste(antigen, "_t_x1",sep="")]), 
         y=x1_elisa[, ag], 
         ylim=c(3,20),
         main=paste("Apac_X1", 
                    microarray_spots$description[microarray_spots$UID==ag], 
                    microarray_spots$INDEX[microarray_spots$UID==ag], sep=", "),
         xlab=paste("log ", antigen, " ELISA titre", sep=""),
         ylab="array intensity")
  }
}
dev.off()

pdf(file="./Figures/Exploratory/0007.KTarray.Apac_X2.scatter_elisa_vs_array.pdf",width=12,height=9)
par(mfcol=c(2,3))
x2_elisa <- elisa[elisa$study=="Apac_X2",]
for(antigen in c("ama", "msp1", "msp2", "csp")){
  microarray_spots <- spots[grep(toupper(antigen), spots$description),]
  for(ag in microarray_spots$UID){
    print(antigen)
    print(ag)
    plot(x=log(x2_elisa[,paste(antigen, "_t_x2",sep="")]), 
         y=x2_elisa[, ag], 
         ylim=c(3,20),
         main=paste("Apac_X2", 
                    microarray_spots$description[microarray_spots$UID==ag], 
                    microarray_spots$INDEX[microarray_spots$UID==ag], sep=", "),
         xlab=paste("log ", antigen, " ELISA titre", sep=""),
         ylab="array intensity")
  }
}
dev.off()

pdf(file="./Figures/Exploratory/0007.KTarray.Apac_X3.scatter_elisa_vs_array.pdf",width=12,height=9)
par(mfcol=c(2,3))
x3_elisa <- elisa[elisa$study=="Apac_X3",]
for(antigen in c("ama", "msp1", "msp2", "csp")){
  microarray_spots <- spots[grep(toupper(antigen), spots$description),]
  for(ag in microarray_spots$UID){
    print(antigen)
    print(ag)
    plot(x=log(x3_elisa[,paste(antigen, "_t_x3",sep="")]), 
         y=x3_elisa[, ag], 
         ylim=c(3,20),
         main=paste("Apac_X3", 
                    microarray_spots$description[microarray_spots$UID==ag], 
                    microarray_spots$INDEX[microarray_spots$UID==ag], sep=", "),
         xlab=paste("log ", antigen, " ELISA titre", sep=""),
         ylab="array intensity")
  }
}
dev.off()
```
