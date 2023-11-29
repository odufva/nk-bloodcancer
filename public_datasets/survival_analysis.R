
# Univariable survival associations for CRISPR screen hits across hematological malignancies

library(RColorBrewer)
library(survival)
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(survMisc)
library(survminer)
library(plyr)
library(dplyr)

# load scripts
source("../CCLE_featurematrix_NK_PRISM/compute.pairwise.R")
source("../CCLE_featurematrix_NK_PRISM/functions_statistics.R")
source("../functions/plotting_functions_may2020.R")
source("../functions/statistics_wrappers_may2020.R")


# load screen results
crispr <- fread("crispr_mageck_combined.txt", data.table = F)

p_threshold <- 0.0001
lfc_threshold_pos <- 0.75
lfc_threshold_neg <- -0.75

activating_lof <- crispr %>%
  filter(!grepl("GOF", cell_line)) %>%
  filter(p < p_threshold, lfc>lfc_threshold_pos) %>%
  select(gene) %>%
  unique() %>%
  tibble::deframe()

activating_gof <- crispr %>%
  filter(grepl("GOF", cell_line)) %>%
  filter(p < p_threshold, lfc<lfc_threshold_neg) %>%
  select(gene) %>%
  unique() %>%
  tibble::deframe()

inhibitory_lof <- crispr %>%
  filter(!grepl("GOF", cell_line)) %>%
  filter(p < p_threshold, lfc<lfc_threshold_neg) %>%
  select(gene) %>%
  unique() %>%
  tibble::deframe()

inhibitory_gof <- crispr %>%
  filter(grepl("GOF", cell_line)) %>%
  filter(p < p_threshold, lfc>lfc_threshold_pos) %>%
  select(gene) %>%
  unique() %>%
  tibble::deframe()

activating <- c(activating_lof, activating_gof)
inhibitory <- c(inhibitory_lof, inhibitory_gof)

genelist  <- unique(c(activating, inhibitory))

#******************************************** Hemap *********************************************

annot = get(load("Hemap_immunology_Annotations.Rdata"))
annot=annot[!is.na(annot$OS_Time),]

# survival time and status
TIME=annot$OS_Time
STATUS=as.numeric(annot$OS_Status)
TIME2=annot$PFS_Time
STATUS2=as.numeric(annot$PFS_Status)

#************************************** Process gene expression data ************************************
data=t(get(load("data9544_with_gene_symbols.RData")))
data=data[,colnames(data)%in%annot$GSM.identifier..sample.]

# data for survival
df=data.frame("time"=annot$OS_Time, "status"=annot$OS_Status, stringsAsFactors = F)

#********************************* just AML *********************************
filterv = annot$subclasses%in%"AML"&annot$CELLS_SORTED==0
logicalv=get.logical(annovector = list(annot$GSE.identifier..experiment.), filterv = filterv)
names(logicalv)=unique(paste(names(logicalv), annot$subclasses[filterv]))

filterv = annot$subclasses%in%"AML"&annot$CELLS_SORTED==0
logicalv2=get.logical(annovector = list(annot$subclasses), filterv = filterv)
names(logicalv2)="Hemap_AML"
logicalv=append(logicalv2, logicalv)
logicalv=logicalv[lapply(logicalv, sum)>2]

# clinical
DATA_clin=data.frame("Age"=annot$AGE, "Age"=annot$AGE, "Gender"=factor(annot$GENDER, levels = c("female", "male")))
aml_res_clin=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_clin,TIME,STATUS, univariate = T, pretty=F))
aml_res_multivar=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_clin,TIME,STATUS, univariate = F, pretty=F))

# genelist
DATA=data.frame(scale(t(data[rownames(data)%in%c(genelist),])))
aml_res_genelist=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# subtype
load("../Bloodcancer_subtypes/subtypes/Hemap_AML_subtypes.Rdata")
coordinates.subtype=coordinates.subtype[!is.na(coordinates.subtype$subtype),]

samples=lapply(unique(coordinates.subtype$subtype), function(type)coordinates.subtype$ID[coordinates.subtype$subtype%in%type])
DATA=data.frame(do.call(cbind, lapply(samples, function(id)annot$GSM.identifier..sample.%in%id))*1)
colnames(DATA)=gsub("-", ".", unique(coordinates.subtype$subtype))

aml_res_subtype=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# make supplementary table, adjusted p-value set here to correct for number of comparisons in total:
result=rbind(aml_res_clin, aml_res_genelist)
result$Adj.P=p.adjust(result$P, method="BH")

result[,2]=prettyNum(result[,2])
result[,3]=prettyNum(result[,3])
result[,4]=prettyNum(result[,4])
result[,5]=prettyNum(result[,5])
result[,6]=prettyNum(result[,6])
result[,8]=prettyNum(result[,8])

# annotate these genes, needed later:
result[,1]=gsub("\\.", "-", result[,1])

genelist_signif=data.frame(result[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%genelist]="CRISPRhit"
genelist_signif$type[genelist_signif[,1]%in%unique(coordinates.subtype$subtype)]="Subtype"
genelist_signif$type[genelist_signif[,1]%in%genelist]="CRISPRhit"
result$Type=genelist_signif$type
result=result[order(result$Type),]

data.table::fwrite(result[,c(1,11,10,2,3,4,5,6,8,9)], "result_Hemap_AML.tsv", sep="\t")
data.table::fwrite(result[result$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "result_Hemap_AML_signif.tsv", sep="\t")

# save cohort and survival
gexp=data.frame(scale(t(data)))

samples=lapply(unique(coordinates.subtype$subtype), function(type)coordinates.subtype$ID[coordinates.subtype$subtype%in%type])
subtypes=do.call(cbind, lapply(samples, function(id)annot$GSM.identifier..sample.%in%id))*1 
colnames(subtypes)=gsub("-", ".", unique(coordinates.subtype$subtype))
samp=data.frame(subtypes, "Age"=annot$AGE, "Gender"=annot$GENDER)

save(list = c("gexp", "TIME", "STATUS", "logicalv", "samp"), file="Hemap_AML_survival_data.Rdata")

#************************************* just MM *************************************
filterv = annot$subclasses%in%"Cancer_Myeloma"&!is.na(annot$MM_ISS)
logicalv=get.logical(annovector = list(annot$GSE.identifier..experiment.), filterv = filterv)
names(logicalv)=unique(paste(names(logicalv), annot$subclasses[filterv]))

# just MM
filterv = annot$subclasses%in%"Cancer_Myeloma"&!is.na(annot$MM_ISS)
logicalv2=get.logical(annovector = list(annot$subclasses), filterv = filterv)
names(logicalv2)="MM_all"
logicalv=append(logicalv2, logicalv)
names(logicalv)=c("Hemap_MM", "GSE19784_Hemap_MM", "GSE16716,GSE24080_Hemap_MM")

data.test=data.frame("time"=TIME[logicalv[[1]]], "status"=STATUS[logicalv[[1]]])

ggsurvplot(survfit(Surv(time, status) ~ 1, data = data.test), 
           xlab = "months", 
           ylab = "Overall survival probability")


DATA=data.frame("ISS3"=(annot$MM_ISS==3)*1, "ISS1"=(annot$MM_ISS==1)*1, "Age"=annot$AGE,
                "Gender"=annot$GENDER, check.names = F)

# clinical
mm_res_clin=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))
mm_res_multivar=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA[,3:5],TIME,STATUS, univariate = F, pretty=F))

# genelist
DATA=data.frame(scale(t(data[rownames(data)%in%genelist,])))
mm_res_genelist=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# subtype
load("Hemap_MM_subtypes.Rdata")

filterv = annot$GSE.identifier..experiment.%in%"GSE16716,GSE24080"&!is.na(annot$MM_ISS)
logicalv=get.logical(annovector = list(annot$GSE.identifier..experiment.), filterv = filterv)
names(logicalv)=c("GSE16716,GSE24080_Hemap_MM")

samples=lapply(unique(coordinates.subtype$subtype), function(type)coordinates.subtype$ID[coordinates.subtype$subtype%in%type])
DATA=data.frame(do.call(cbind, lapply(samples, function(id)annot$GSM.identifier..sample.%in%id))*1)
colnames(DATA)=gsub("-", "_", unique(coordinates.subtype$subtype))

mm_res_subtype=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# make supplementary table, adjusted p-value set here to correct for number of comparisons in total:
result=rbind(mm_res_clin, mm_res_genelist, mm_res_subtype)
result$Adj.P=p.adjust(result$P, method="BH")

result[,2]=prettyNum(result[,2])
result[,3]=prettyNum(result[,3])
result[,4]=prettyNum(result[,4])
result[,5]=prettyNum(result[,5])
result[,6]=prettyNum(result[,6])
result[,8]=prettyNum(result[,8])

# annotate these genes, needed later:
genelist_signif=data.frame(result[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%gsub("-", "_", unique(coordinates.subtype$subtype))]="Subtype"
genelist_signif$type[genelist_signif[,1]%in%gsub("\\-", ".", genelist)]="CRISPRhit"
result$Type=genelist_signif$type
result=result[order(result$Type),]
result[,1]=gsub("\\.", "-", result[,1])

data.table::fwrite(result[,c(1,11,10,2,3,4,5,6,8,9)], "result_Hemap_MM.tsv", sep="\t")
data.table::fwrite(result[result$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "result_Hemap_MM_signif.tsv", sep="\t")


# save cohort and survival
gexp=data.frame(scale(t(data)))

samples=lapply(unique(coordinates.subtype$subtype), function(type)coordinates.subtype$ID[coordinates.subtype$subtype%in%type])
subtype=do.call(cbind, lapply(samples, function(id)annot$GSM.identifier..sample.%in%id))*1
colnames(subtype)=c(gsub("-", ".", unique(coordinates.subtype$subtype)), colnames(subtype))
samp=data.frame(subtype, "Age"=annot$AGE, "Gender"=annot$GENDER, "ISS3"=(annot$MM_ISS==3)*1, "ISS1"=(annot$MM_ISS==1)*1)

save(list = c("gexp", "TIME", "STATUS", "logicalv", "samp"), file="Hemap_MM_survival_data.Rdata")

#***************************** just DLBCL ************************************

filterv = annot$subclasses%in%"BCL_DLBCL"&!is.na(annot$dlbcl_ipi)
logicalv=get.logical(annovector = list(annot$GSE.identifier..experiment.), filterv = filterv)
names(logicalv)=unique(paste(names(logicalv), annot$subclasses[filterv]))

filterv = annot$subclasses%in%"BCL_DLBCL"
logicalv2=get.logical(annovector = list(annot$subclasses), filterv = filterv)
names(logicalv2)="BCL_DLBCL_all"

logicalv=append(logicalv2, logicalv)
logicalv=logicalv[lapply(logicalv, sum)>2]

RCHOP=list(annot$Chemotherapy_RCHOP==1&!is.na(STATUS)&!TIME==0)
CHOP=list(annot$Chemotherapy_CHOP==1&!is.na(STATUS)&!TIME==0)
names(RCHOP)="Hemap_DLBCL_RCHOP"
names(CHOP)="Hemap_DLBCL_CHOP"

logicalv <- append(logicalv, append(RCHOP, CHOP))

# clinical
DATA=data.frame("IPI_0to1"=(annot$dlbcl_ipi%in%c(0,1))*1,
                "IPI_4to5"=(annot$dlbcl_ipi%in%c(4,5))*1,
                "ABC"=(grepl("ABC", annot$tbLY))*1,
                "GCB"=(grepl("GCB", annot$tbLY))*1,check.names = F)

dlbcl_res_clin=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# genelist
DATA=data.frame(scale(t(data[rownames(data)%in%c(genelist),])))
dlbcl_res_genelist=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))


# make supplementary table, adjusted p-value set here to correct for number of comparisons in total:
result=rbind(dlbcl_res_clin, dlbcl_res_genelist)
result=result[result$Name%in%c("Hemap_DLBCL_RCHOP"),]
result$Adj.P=p.adjust(result$P, method="BH")

result[,2]=prettyNum(result[,2])
result[,3]=prettyNum(result[,3])
result[,4]=prettyNum(result[,4])
result[,5]=prettyNum(result[,5])
result[,6]=prettyNum(result[,6])
result[,8]=prettyNum(result[,8])

# annotate these genes, needed later:
genelist_signif=data.frame(result[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%genelist]="CRISPRhit"
genelist_signif$type[result$Feature%in%c("ABC", "GCB")]="Subtype"

result$Type=genelist_signif$type
result=result[order(result$Type),]
result[,1]=gsub("\\.", "-", result[,1])

data.table::fwrite(result[,c(1,11,10,2,3,4,5,6,8,9)], "result_Hemap_DLBCL_RCHOP.tsv", sep="\t")
data.table::fwrite(result[result$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "result_Hemap_DLBCL_RCHOP_signif.tsv", sep="\t")


# make supplementary table, adjusted p-value set here to correct for number of comparisons in total:
result=rbind(dlbcl_res_clin, dlbcl_res_genelist)
result=result[result$Name%in%c("Hemap_DLBCL_CHOP"),]
result$Adj.P=p.adjust(result$P, method="BH")

result[,2]=prettyNum(result[,2])
result[,3]=prettyNum(result[,3])
result[,4]=prettyNum(result[,4])
result[,5]=prettyNum(result[,5])
result[,6]=prettyNum(result[,6])
result[,8]=prettyNum(result[,8])

# annotate these genes, needed later:
genelist_signif=data.frame(result[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%genelist]="CRISPRhit"
genelist_signif$type[result$Feature%in%c("ABC", "GCB")]="Subtype"

result$Type=genelist_signif$type
result=result[order(result$Type),]
result[,1]=gsub("\\.", "-", result[,1])

data.table::fwrite(result[,c(1,11,10,2,3,4,5,6,8,9)], "result_Hemap_DLBCL_CHOP.tsv", sep="\t")
data.table::fwrite(result[result$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "result_Hemap_DLBCL_CHOP_signif.tsv", sep="\t")

# save cohort and survival
gexp=data.frame(scale(t(data)))

samp=data.frame("Age"=annot$AGE, "Gender"=annot$GENDER, "IPI_0to1"=(annot$dlbcl_ipi%in%c(0,1))*1, "IPI_4to5"=(annot$dlbcl_ipi%in%c(4,5))*1,"ABC"=(grepl("ABC", annot$tbLY))*1, "GCB"=(grepl("GCB", annot$tbLY))*1)

save(list = c("gexp","samp", "TIME", "STATUS", "logicalv"), file="Hemap_DLBCL_survival_data.Rdata")


#********************************************************* Compass MM *********************************************************
fm=get(load("MM_COMPASS_FM.Rdata"))
annot=get(load("MM_COMPASS_ANNOT.Rdata"))

annot=annot[match(colnames(fm), rownames(annot)),]

data=fm[grepl("N:GEXP:", rownames(fm)),]
rownames(data)=gsub("N:GEXP:", "", rownames(data))

data.mut=fm[grepl("B:GNAB:", rownames(fm)),]
rownames(data.mut)=gsub("B:GNAB:", "", rownames(data.mut))
data.mut.filt=data.mut[!rowSums(data.mut, na.rm = T)<10,]

TIME=as.numeric(fm["N:CLIN:OS",])
STATUS=as.numeric(fm["B:CLIN:STATUS",])

STATUS[TIME>1825&!is.na(STATUS)&STATUS==1]=0 # change to 4year survival, 5 year sharp drop
TIME[TIME>1825]=1825

# plot overall survival
r=lapply(seq(logicalv), function(i){
  data.test=data.frame("time"=TIME[logicalv[[i]]], "status"=STATUS[logicalv[[i]]])
  
  ggsurvplot(survfit(Surv(time, status) ~ 1, data = data.test), 
             xlab = "months", 
             ylab = "Overall survival probability",title=(names(logicalv)[i])
  )
})

names(r)=names(logicalv)

pdf("CoMMpass_MM_cohorts.pdf", width =5, height = ceiling(length(r)/2)*2.75)
plots.together=arrange_ggsurvplots(r, print = TRUE, ncol = 2, nrow = ceiling(length(r)/2))
dev.off()

# univariable analysis for clinical

logicalv=list(!is.na(data["CFLAR",]))
names(logicalv)="CoMMpass"

bortezomib=list(annot$D_PT_therclass=="Bortezomib-based")
bortezomib_IMID=list(annot$D_PT_therclass=="combined bortezomib/IMIDs-based")
names(bortezomib)="CoMMpass_bortezomib"
names(bortezomib_IMID)="CoMMpass_bortezomib_IMID"

logicalv <- append(logicalv, append(bortezomib, bortezomib_IMID))

DATA=data.frame("ISS1"=as.numeric(fm["B:CLIN:R_ISS_1",]),
                "ISS2"=as.numeric(fm["B:CLIN:R_ISS_2",]),
                "ISS3"=as.numeric(fm["B:CLIN:R_ISS_3",]), stringsAsFactors = F)

commpass_res=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# univariable analysis for genelist

genelist = genelist[genelist %in% rownames(data)]


DATA=data.frame(scale(t(data[rownames(data)%in%c(genelist),])))
commpass_genelist=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA, TIME, STATUS, univariate = T, pretty = F))


# subtype:
load("CoMMpass_MM_subtypes.Rdata")
coordinates.subtype=coordinates.subtype[match(colnames(fm), coordinates.subtype$ID),]

DATA_subtypes=data.frame(do.call(cbind, get.logical(list(coordinates.subtype$subtype)))*1)
commpass_subtype=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_subtypes,TIME,STATUS, univariate = T, pretty=F))

# make supplementary table, adjusted p-value set here to correct for number of comparisons in total:
result=rbind(commpass_res, commpass_genelist, commpass_subtype)
result$Adj.P=p.adjust(result$P, method="BH")

# annotate these genes, needed later:
result[,2]=prettyNum(result[,2])
result[,3]=prettyNum(result[,3])
result[,4]=prettyNum(result[,4])
result[,5]=prettyNum(result[,5])
result[,6]=prettyNum(result[,6])
result[,8]=prettyNum(result[,8])

# annotate these genes, needed later:
genelist_signif=data.frame(result[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%genelist]="CRISPRhit"
genelist_signif$type[genelist_signif[,1]%in%commpass_subtype$Feature]="Subtype"

result$Type=genelist_signif$type
result=result[order(result$Type),]
result[,1]=gsub("\\.", "-", result[,1])

data.table::fwrite(result[,c(1,11,10,2,3,4,5,6,8,9)], "result_CoMMpass.tsv", sep="\t")
data.table::fwrite(result[result$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "result_CoMMpass_signif.tsv", sep="\t")

# save cohort and survival
gexp=data.frame(scale(t(data)))

samp=cbind(data.frame("ISS1"=as.numeric(fm["B:CLIN:R_ISS_1",]),"ISS2"=as.numeric(fm["B:CLIN:R_ISS_2",]),"ISS3"=as.numeric(fm["B:CLIN:R_ISS_3",]), stringsAsFactors = F), DATA_subtypes)

save(list = c("gexp", "samp", "TIME", "STATUS", "logicalv"), file="CoMMpass_survival_data.Rdata")


#********************************************************* DLBCL GSE985888 *********************************************************
fm=get(load("../Chapuy_DLBCL_NK_PRISM/GSE98588_fm.Rdata"))
annot=get(load("GSE98588_annot.Rdata"))

annot=annot[match(colnames(fm), annot$individual_id),]

data=fm[grepl("N:GEXP:", rownames(fm)),]
rownames(data)=gsub("N:GEXP:", "", rownames(data))

data.mut=fm[grepl("B:GNAB:", rownames(fm)),]
rownames(data.mut)=gsub("B:GNAB:", "", rownames(data.mut))
data.mut.filt=data.mut[!rowSums(data.mut, na.rm = T)<10,]

TIME=as.numeric(fm["N:CLIN:OS",])
STATUS=as.numeric(fm["B:CLIN:OS_STAT",])


# univariable analysis for clinical

logicalv=list(!is.na(data["CFLAR",]))
names(logicalv)="DLBCL_GSE98588"

DATA=data.frame("IPI_0to1"=(annot$IPI%in%c(0,1))*1,
                "IPI_4to5"=(annot$IPI%in%c(4,5))*1,
                "ABC"=(grepl("ABC", annot$COO_byGEP))*1,
                "GCB"=(grepl("GCB", annot$COO_byGEP))*1,check.names = F)

chapuy_res=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# univariable analysis for genelist

genelist = genelist[genelist %in% rownames(data)]

#logicalv=list(!is.na(data["CFLAR",]) & annot$subtype.cluster %in% c("WHSC1_FGFR3_Ig"))

DATA=data.frame(scale(t(data[rownames(data)%in%c(genelist),])))
chapuy_genelist=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA, TIME, STATUS, univariate = T, pretty = F))


# make supplementary table, adjusted p-value set here to correct for number of comparisons in total:
result=rbind(chapuy_res, chapuy_genelist)
result$Adj.P=p.adjust(result$P, method="BH")

# annotate these genes, needed later:
result[,2]=prettyNum(result[,2])
result[,3]=prettyNum(result[,3])
result[,4]=prettyNum(result[,4])
result[,5]=prettyNum(result[,5])
result[,6]=prettyNum(result[,6])
result[,8]=prettyNum(result[,8])

# annotate these genes, needed later:
genelist_signif=data.frame(result[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%genelist]="CRISPRhit"

result$Type=genelist_signif$type
result=result[order(result$Type),]
result[,1]=gsub("\\.", "-", result[,1])

data.table::fwrite(result[,c(1,11,10,2,3,4,5,6,8,9)], "result_GSE98588.tsv", sep="\t")
data.table::fwrite(result[result$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "result_GSE98588_signif.tsv", sep="\t")

# save cohort and survival
gexp=data.frame(scale(t(data)))

samp=data.frame("Age"=annot$Age.at.first.diagnosis, "Gender"=annot$Gender, "IPI_0to1"=(annot$IPI%in%c(0,1))*1, "IPI_4to5"=(annot$IPI%in%c(4,5))*1,"ABC"=(grepl("ABC", annot$COO_byGEP))*1, "GCB"=(grepl("GCB", annot$COO_byGEP))*1)

save(list = c("gexp", "samp", "TIME", "STATUS", "logicalv"), file="GSE98588_survival_data.Rdata")



#********************************************************* Reddy DLBCL *********************************************************
fm=get(load("../Reddy_DLBCL/REDDY_DLBCL_fm.Rdata"))
annot=get(load("REDDY_DLBCL_annot.Rdata"))

data=data.matrix(fm[grepl("GEXP", rownames(fm)),])
rownames(data)=gsub("N:GEXP:", "", rownames(data))

TIME=as.numeric(annot$Overall.Survival.years)
STATUS=(as.numeric(annot$Censored)==0)*1 # 0 indicates no censoring, meaning that the death was observed, whereas a 1 indicates that the patient was alive

l.regulon.gene=regulon.feats(fm, genelist)

# all costim-cytolytic-CGA correlated mutations:
genelist_feats=data.matrix(fm[rownames(fm)%in%unlist(l.regulon.gene),])
rownames(genelist_feats)=sapply(rownames(genelist_feats), function(n){
  if(!grepl("CNVR", n))return(n)
  a=names(l.regulon.gene)[sapply(l.regulon.gene, function(a)any(a%in%n))]
  # if(length(a))paste0(n,"@", paste(a, collapse=","))
})

DATAmut=data.frame(t(genelist_feats[grepl("GNAB|CNVR", rownames(genelist_feats)),]))

logicalv=list("Reddy_DLBCL"=!is.na(data[1,]))

Reddy_DLBCL_res=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))
Reddy_DLBCL_multivar=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = F, pretty=F))

# Clinical
DATA_clin=data.frame("IPI_0to1"=annot$IPI%in%c(0,1)*1,"IPI_4to5"=annot$IPI%in%c(4,5)*1, "ABC"=(annot$ABC.GCB..RNAseq.=="ABC")*1 ,"GCB"=(annot$ABC.GCB..RNAseq.=="GCB")*1,  stringsAsFactors = F)
Reddy_DLBCL_clin=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_clin,TIME,STATUS, univariate = T, pretty=F))

# immuno-editing:
Reddy_DLBCL_immunoediting=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATAmut,TIME,STATUS, univariate = T, pretty=F))
Reddy_DLBCL_immunoediting$Feature=gsub("B.GNAB.", "MUT:", Reddy_DLBCL_immunoediting$Feature)
Reddy_DLBCL_immunoediting$Feature=gsub("HLA.", "HLA-", Reddy_DLBCL_immunoediting$Feature)
Reddy_DLBCL_immunoediting$Feature=gsub("N.CNVR.", "", Reddy_DLBCL_immunoediting$Feature)

# genelist
DATA=data.frame(scale(t(data[rownames(data)%in%c(genelist),])))
Reddy_DLBCL_genelist=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# make supplementary table, adjusted p-value set here to correct for number of comparisons in total:
result=rbind(Reddy_DLBCL_res,Reddy_DLBCL_clin[!is.na(Reddy_DLBCL_clin$`exp(coef)`),], Reddy_DLBCL_genelist, Reddy_DLBCL_immunoediting) 
result=result[result$Name%in%"Reddy_DLBCL",]
result$Adj.P=p.adjust(result$P, method="BH")


result[,2]=prettyNum(result[,2])
result[,3]=prettyNum(result[,3])
result[,4]=prettyNum(result[,4])
result[,5]=prettyNum(result[,5])
result[,6]=prettyNum(result[,6])
result[,8]=prettyNum(result[,8])

# annotate these genes, needed later:
genelist_signif=data.frame(result[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%genelist]="CRISPRhit"
genelist_signif$type[genelist_signif[,1]%in%Reddy_DLBCL_immunoediting$Feature]="Immune Editing mutation"
genelist_signif$type[result$Feature%in%c("ABC", "GCB")]="Subtype"

result$Type=genelist_signif$type
result=result[order(result$Type),]
result[,1]=gsub("\\.", "-", result[,1])

data.table::fwrite(result[,c(1,11,10,2,3,4,5,6,8,9)], "result_Reddy_DLBCL.tsv", sep="\t")
data.table::fwrite(result[result$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "result_Reddy_DLBCL_signif.tsv", sep="\t")

# save cohort and survival
gexp=data.frame(scale(t(data)))

samp=DATA_clin

save(list = c("gexp", "samp", "TIME", "STATUS", "logicalv"), file="Reddy_DLBCL_survival_data.Rdata")


#************************************ beatAML:
load("BeatAML_fm.Rdata")
annot=get(load("BeatAML_fm_annot.Rdata"))

gexp=fm[grepl("GEXP", rownames(fm)),]
rownames(gexp)=gsub("N:GEXP:", "", rownames(gexp))

annot$vitalStatus[annot$vitalStatus=="Unknown"]=NA
TIME=as.numeric(annot$overallSurvival)
STATUS=as.numeric(annot$vitalStatus=="Dead")

STATUS[TIME>1825&!is.na(STATUS)&STATUS==1]=0 # change to 5year survival
TIME[TIME>1825]=1825
TIME[TIME==0]=0.1

filtv=annot$specimenType%in%"Bone Marrow Aspirate"&!(is.na(TIME)|is.na(STATUS)|is.na(annot$TCGA_coord))&annot$causeOfDeath%in%c("Alive", "Dead-Disease", "Dead-Unknown")
filtv=annot$specimenType%in%"Bone Marrow Aspirate"&!(is.na(annot$TCGA_coord))&annot$causeOfDeath%in%c("Alive", "Dead-Disease", "Dead-Unknown")
filtv2=annot$specimenType%in%"Bone Marrow Aspirate"&!(is.na(annot$TCGA_coord))&annot$TCGA_coord%in%c("CMP-like","MDS-like","Monocyte-like")&annot$causeOfDeath%in%c("Alive", "Dead-Disease", "Dead-Unknown")
filtv3=annot$specimenType%in%"Bone Marrow Aspirate"&!(is.na(annot$TCGA_coord))&annot$TCGA_coord%in%c("MDS-like")&annot$causeOfDeath%in%c("Alive", "Dead-Disease", "Dead-Unknown")
filtv4=annot$specimenType%in%"Bone Marrow Aspirate"&!(is.na(annot$TCGA_coord))&annot$TCGA_coord%in%c("CMP-like")&annot$causeOfDeath%in%c("Alive", "Dead-Disease", "Dead-Unknown")
filtv5=annot$specimenType%in%"Bone Marrow Aspirate"&!(is.na(annot$TCGA_coord))&annot$TCGA_coord%in%c("Monocyte-like")&annot$causeOfDeath%in%c("Alive", "Dead-Disease", "Dead-Unknown")

logicalv=list("BeatAML"=!(is.na(TIME)|is.na(STATUS)|is.na(annot$TCGA_coord))&annot$causeOfDeath%in%c("Alive", "Dead-Disease", "Dead-Unknown"), "beatAML_BMonly"=filtv)#, "normal_karyotype"=filtv2, "MDS-like"=filtv3,  "CMP-like"=filtv4,  "Monocyte-like"=filtv5)

data.test=data.frame("time"=TIME[logicalv[[1]]]*0.0328767, "status"=STATUS[logicalv[[1]]])

ggsurvplot(survfit(Surv(time, status) ~ 1, data = data.test), 
           xlab = "months", 
           ylab = "Overall survival probability")

# CRISPR hits

genelist = genelist[genelist %in% rownames(gexp)]

DATA=data.frame(scale(t(gexp[rownames(gexp)%in%c(genelist),])))
beatAML_genelist=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# Clinical:
DATA_clin=data.frame(t(fm[grepl("ELN2017", rownames(fm)),]),"Age"=as.numeric(fm[grepl("ageAtDiagnosis", rownames(fm)),,drop=F]), t(fm[grepl("BM_Transplant", rownames(fm)),,drop=F]),t(fm[grepl("B:CLIN:priorMDS_TRUE", rownames(fm)),,drop=F]), t(fm[grepl("in_PB|in_BM|B:CLIN:is|finalFusion_Complex", rownames(fm)),,drop=F]))
DATA_clin$Age[is.na(DATA_clin$Age)]=median(DATA_clin$Age, na.rm = T) #one observation --> set to median
DATA_clin=DATA_clin[,colSums(DATA_clin, na.rm = T)>3]

colnames(DATA_clin)=gsub("..CLIN.", "", colnames(DATA_clin))
beatAML_clin=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_clin,TIME,STATUS, univariate = T, pretty=F))

# subtype:
load("BeatAML_subtypes.Rdata")
coordinates.subtype$subtype[coordinates.subtype$subtype%in%"Progenitor-like"]="CMP-like"
coordinates.subtype=coordinates.subtype[match(colnames(fm), coordinates.subtype$ID),]

DATA_subtypes=data.frame(do.call(cbind, get.logical(list(coordinates.subtype$subtype)))*1)

beataml_res_subtype=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_subtypes,TIME,STATUS, univariate = T, pretty=F))

fun.kapplanMeier(TIME[logicalv[[1]]], STATUS[logicalv[[1]]],GROUPS=coordinates.subtype$subtype[logicalv[[1]]], MONTHS=T, PVAL=1, INDIVIDUAL_GROUPS=F,LWD = 1, NAME = "Prognostic Index - validation")

# make supplementary table, adjusted p-value set here to correct for number of comparisons in total:
result=rbind(beatAML_genelist,beatAML_clin,beataml_res_subtype) 
result=result[result$Name%in%"BeatAML",]
result$Adj.P=p.adjust(result$P, method="BH")

result[,2]=prettyNum(result[,2])
result[,3]=prettyNum(result[,3])
result[,4]=prettyNum(result[,4])
result[,5]=prettyNum(result[,5])
result[,6]=prettyNum(result[,6])
result[,8]=prettyNum(result[,8])

# annotate these genes, needed later:
genelist_signif=data.frame(result[,1], type="", stringsAsFactors = F)
genelist_signif=data.frame(result[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%genelist]="CRISPRhit"
genelist_signif$type[genelist_signif[,1]%in%gsub("\\-", ".", unique(coordinates.subtype$subtype))]="Subtype"

result$Type=genelist_signif$type
result=result[order(result$Type),]
result[,1]=gsub("\\.", "-", result[,1])

data.table::fwrite(result[,c(1,11,10,2,3,4,5,6,8,9)], "result_beatAML.tsv", sep="\t")
data.table::fwrite(result[result$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "result_beatAML_signif.tsv", sep="\t")

# save cohort and survival
gexp=data.frame(scale(t(gexp)))

samp=cbind(DATA_clin, DATA_subtypes)

save(list = c("gexp", "samp", "TIME", "STATUS", "logicalv"), file="BeatAML_survival_data.Rdata")


#******************************************** TCGA AML *****************************************************
fm_org=get(load("TCGA_AML_FM_DUFVA.Rdata"))
fm=fm_org[,!is.na(fm_org["N:SAMP:CytolyticScore",])]

risks=c("N:CLIN:Age:::::",
        "C:CLIN:acute_myeloid_leukemia_calgb_cytogenetics_risk_category:::::" ,
        "C:CLIN:FISH_test_component:::::",
        "B:GNAB:NPM1:chr5:170814708:170837888:+:y_n_somatic",
        "B:GNAB:FLT3:chr13:28577411:28674729:-:y_n_somatic",
        "B:GNAB:CEBPA:chr19:33790840:33793430:-:y_n_somatic",
        "B:GNAB:TP53:chr17:7565097:7590863:-:y_n_somatic")

df=t(fm[rownames(fm)%in%risks,])
colnames(df)=do.call(rbind, strsplit(colnames(df), ":"))[,3]

data=data.matrix(fm[grepl("GEXP", rownames(fm)),])
rownames(data)=make.unique(do.call(rbind, strsplit(rownames(data), ":"))[,3])


OS=as.numeric(fm["N:CLIN:OS.months..3.31.12:::::",])
TIME=OS
STATUS=as.numeric(fm["C:CLIN:vital_status_TCGA_paper:::::",]=="DECEASED")
PFS=as.numeric(fm["N:CLIN:EFS.months....4.30.13:::::",])


logicalv=list("TCGA_AML"=rep(T, dim(DATA)[1]))

# CRISPR hits

genelist = genelist[genelist %in% rownames(data)]

DATA=data.frame(scale(t(data[rownames(data)%in%c(genelist),])))
TCGA_AML_genelist=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))


# clinical
DATA_clin=data.frame("Age"=scale(as.numeric(fm["N:CLIN:Age:::::",])), "Blast.percentage"=scale(as.numeric(fm["N:CLIN:X.BM.Blast:::::",])), "Gender"=as.character(fm["C:CLIN:Sex:::::",]))
TCGA_AML_clin=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_clin,TIME,STATUS, univariate = T, pretty=F))

# subtype:
load("TCGA_AML_subtypes.Rdata")

DATA_subtypes=data.frame(do.call(cbind, get.logical(list(coordinates.subtype$subtype)))*1)

TCGA_aml_res_subtype=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_subtypes,TIME,STATUS, univariate = T, pretty=F))

# make supplementary table, adjusted p-value set here to correct for number of comparisons in total:
result=rbind(TCGA_AML_genelist, TCGA_aml_res_subtype, TCGA_AML_clin)
result$Adj.P=p.adjust(result$P, method="BH")

result[,2]=prettyNum(result[,2])
result[,3]=prettyNum(result[,3])
result[,4]=prettyNum(result[,4])
result[,5]=prettyNum(result[,5])
result[,6]=prettyNum(result[,6])
result[,8]=prettyNum(result[,8])

# annotate these genes, needed later:
result[,1]=gsub("\\.", "-", result[,1])

genelist_signif=data.frame(result[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%genelist]="CRISPRhit"
genelist_signif$type[genelist_signif[,1]%in%unique(coordinates.subtype$subtype)]="Subtype"

result$Type=genelist_signif$type
result=result[order(result$Type),]

data.table::fwrite(result[,c(1,11,10,2,3,4,5,6,8,9)], "result_TCGA_AML.tsv", sep="\t")
data.table::fwrite(result[result$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "result_TCGA_AML_signif.tsv", sep="\t")

# save cohort and survival
gexp=data.frame(scale(t(data)))

samp=cbind(DATA_clin, DATA_subtypes)

save(list = c("gexp", "samp", "TIME", "STATUS", "logicalv"), file="TCGA_AML_survival_data.Rdata")




