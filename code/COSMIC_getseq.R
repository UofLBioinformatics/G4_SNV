#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)

library(dplyr)
library(purrr)
library(ggplot2)

require(Biostrings) 
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)


mcsapply<-function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, mc.preschedule = TRUE,
 mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L), 
  mc.cleanup = TRUE, mc.allow.recursive = TRUE, affinity.list = NULL )
{
  answer <- mclapply(X = X, FUN = FUN, ...,mc.preschedule = mc.preschedule, 
mc.set.seed = mc.set.seed, mc.silent = mc.silent, mc.cores = mc.cores, 
  mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive, affinity.list = affinity.list)
  if (USE.NAMES && is.character(X) && is.null(names(answer))) 
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer)) 
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}






getseq=function(filename,typeofvar){

SubDir <- paste0(mainDir,added_name_for_output_file,"/",typeofvar,"/")
#print(SubDir)


dir.create(file.path( SubDir ), showWarnings = FALSE,recursive=TRUE)


setwd(file.path( SubDir ))
start_time <- Sys.time()

genome<-BSgenome.Hsapiens.UCSC.hg38

VCF<-fread(file=filename)

colnames(VCF)[1]<-"chr"
chromosome<-c(1:22,"X","Y")
VCF1<-VCF[VCF$chr %in% chromosome,]
VCF1$chr<-paste0("chr",VCF1$chr)


g4_maf<-VCF1
g4_maf$start<-g4_maf$POS-30

g4_maf$end<-g4_maf$POS+30
g4_maf$loc_start<-g4_maf$POS
g4_maf$loc_end<-g4_maf$POS

g4_maf$snp_left<-abs(as.numeric(g4_maf$start)-as.numeric(g4_maf$loc_start))
g4_maf$snp_right<-as.numeric(g4_maf$end)-as.numeric(g4_maf$loc_end)

g4_maf<-g4_maf[((nchar(g4_maf$REF)==1) | (nchar(g4_maf$ALT)==1)),]
seq_g4_1 <- BSgenome::getSeq(genome, names = g4_maf$chr,
                               start = (as.numeric(g4_maf$start)+1),
                               end = (as.numeric(g4_maf$start)+as.numeric(g4_maf$snp_left)-1))

seq_g4_2 <- getSeq(genome, names = g4_maf$chr,
                               start = (as.numeric(g4_maf$loc_start)+1),
                               end = (as.numeric(g4_maf$loc_start)+as.numeric(g4_maf$snp_right)+1))
  
  seq_g4_full <- getSeq(genome, names = g4_maf$chr,
                               start = (as.numeric(g4_maf$start)+1),
                               end = (as.numeric(g4_maf$end)+1))
  




head(as.data.frame(seq_g4_1))
head(as.data.frame(seq_g4_2))
head(as.data.frame(seq_g4_full))
g4_seq<-(seq_g4_full)
g4_maf <- g4_maf %>% mutate(G4_id = row_number())



g4_maf_s<-g4_maf
seq_left<-as.data.frame(seq_g4_1)
colnames(seq_left)<-"leftseq"
#last letter of the left seq is the position in question
seq_right<-as.data.frame(seq_g4_2)



rm(seq_g4_1,seq_g4_2)
colnames(seq_right)<-"rightseq"
seq_g4_full<-as.data.frame(seq_g4_full)
colnames(seq_g4_full)<-"seq_g4_full"

seq<-cbind(as.data.frame(seq_left),as.data.frame(g4_maf_s),as.data.frame(seq_right),as.data.frame(seq_g4_full))


seq$G4_ref<-paste0(seq$leftseq,seq$REF,seq$rightseq)
seq$G4_alt<-paste0(seq$leftseq,seq$ALT,seq$rightseq)

all_1<-seq

identical(seq$G4_ref,seq$seq_g4_full)


i=1
seq$name<-paste0(seq$chr,"_",seq$POS,"_","G4_id")
seq.list <- seq[,c("G4_id","name","G4_ref","G4_alt","seq_g4_full","REF","ALT","loc_start")]

out<-paste0(SubDir,typeofvar,".rds")
saveRDS(seq.list,out)


}





filename1<-"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/NCV_indelfiltered.recode.vcf"
typeofvar_1<-"NCV"

filename2<-"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/coding_indelfiltered.recode.vcf"
typeofvar_2<-"coding"

filename3<-"/bio/home/goonerrn/g_quadruplex/annotation_database/SNP/ClinVar/clinvar_SNV.vcf"
typeofvar_2<-"clinvar"

filename4<-"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/1000genomephase3.g4.vcf/1000genomephase3_combined.bed"
VCF<-fread(file=filename4)
VCF<-VCF[,1:9]
colnames(VCF)<-c("chr","POS","V3","REF","ALT","V6","V7","V8","V9")
fwrite(VCF,paste0(filename4,"_edited"),sep="\t",quote=F,row.names=F)
filename4<-"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/1000genomephase3.g4.vcf/1000genomephase3_combined.bed_edited"
#after this editing

typeofvar_2<-"vcf1000genome"
basena<-basename(filename4)

basena<-gsub(pattern = "\\.vcf$", "", basena)


added_name_for_output_file<-args[1]



mainDir<-"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/"
SubDir0 <- paste0(mainDir,added_name_for_output_file,"/")
dir.create(file.path( SubDir0 ), showWarnings = FALSE,recursive=TRUE)


system.time(getseq(filename1,typeofvar_1))
system.time(getseq(filename2,typeofvar_2))
system.time(getseq(filename3,typeofvar_2))

system.time(getseq(filename4,typeofvar_2))









#taken from G4hunter paper
#source("function_G4Hunter.r")
##g4hunter.input1.2<-split(seq.list,1:NROW(seq.list))
#g4hunter.input1.2<-seq.list
#mclapply(g4hunter.input1.2, 
#for (i in 1:nrow(g4hunter.input1.2)){
#  g4hunter.input1.2[i,9]<-(sapply(g4hunter.input1.2[i,3],function(x) signif(G4Hscore(x),3)))%>% unlist()
#  #g4hunter.input1.2$<-cbind(g4hunter.input1.2,sapply(g4hunter.input1.2[i,5],function(x) signif(G4Hscore(x),3)))
#  g4hunter.input1.2[i,10] <- sapply(g4hunter.input1.2[i,3],QPtest)%>% unlist()
#  g4hunter.input1.2[i,11]<- sapply(g4hunter.input1.2[i,3],QPtest,lmax=7)%>% unlist()
#  g4hunter.input1.2[i,12] <-sapply(g4hunter.input1.2[i,3],QPtest,gmin=3,lmax=7)%>% unlist()
#  
#  g4hunter.input1.2[i,13]<-sapply(g4hunter.input1.2[i,4],function(x) signif(G4Hscore(x),3)) %>% unlist()
#  g4hunter.input1.2[i,14] <-  sapply( g4hunter.input1.2[i,4],QPtest)%>% unlist()
#  g4hunter.input1.2[i,15] <- sapply(g4hunter.input1.2[i,4],QPtest,lmax=7)%>% unlist()
#  g4hunter.input1.2[i,16] <- sapply(g4hunter.input1.2[i,4],QPtest,gmin=3,lmax=7)%>% unlist()
#

#}


#end_time <- Sys.time()
#print(difftime(end_time, start_time, units = "secs")[[1]])


dayloop2 <- function(temp){
  temp <- lapply(1:nrow(temp), function(i) {
    cat(round(i/nrow(temp)*100,2),"%    \r")
    #do stuff
  })
  return(temp)
}
