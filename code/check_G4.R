

#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(parallel)
library(dplyr)
library(purrr)
library(ggplot2)

require(Biostrings) 
#library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
start_time<-Sys.time()

args <- commandArgs(trailingOnly = TRUE)
numrow<-args[1]
print(numrow)
numrow<-as.numeric(numrow)
codingseq<-readRDS("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/coding/coding_indelfiltered.recode.rds")
#clinvarseq<-readRDS("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/clinvar.rds")
source("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/function_G4Hunter.r")
print("read sequence")
#codingseq<-clinvarseq
codingseq=codingseq[sample(nrow(codingseq),replace=F,size=numrow*nrow(codingseq)),]


####g4hunter.input1.2<-split(seq.list,1:NROW(seq.list))
start_time<-Sys.time()

#print(head(codingseq))


#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
dog4hunter<-function(x,out,path){
seq.list1<-x
g4hunter.input1.2<-setDT(seq.list1)


V9<-   g4hunter.input1.2[,signif(G4Hscore(G4_ref),3), by = seq_len(nrow(g4hunter.input1.2))]%>% select(V1)
    
V10<- g4hunter.input1.2[,QPtest(G4_ref), by = seq_len(nrow(g4hunter.input1.2))]%>% select(V1)

V11<- g4hunter.input1.2[, QPtest(G4_ref,lmax=7), by = seq_len(nrow(g4hunter.input1.2))]%>% select(V1)

V12<-  g4hunter.input1.2[, QPtest(G4_ref,gmin=3,lmax=7), by = seq_len(nrow(g4hunter.input1.2))]%>% select(V1)

V13<-   g4hunter.input1.2[,signif(G4Hscore(G4_alt),3), by = seq_len(nrow(g4hunter.input1.2))]%>% select(V1)
    
V14<- g4hunter.input1.2[,QPtest(G4_alt), by = seq_len(nrow(g4hunter.input1.2))]%>% select(V1)
V15<- g4hunter.input1.2[, QPtest(G4_alt,lmax=7), by = seq_len(nrow(g4hunter.input1.2))]%>% select(V1)
V16<-  g4hunter.input1.2[, QPtest(G4_alt,gmin=3,lmax=7), by = seq_len(nrow(g4hunter.input1.2))]%>% select(V1)


V1<-g4hunter.input1.2[,1]
V2<-g4hunter.input1.2[,2]
V3<-g4hunter.input1.2[,3]
V4<-g4hunter.input1.2[,4]

result<-c(V1,V2,V3,V4,V9,V10,V11,V12,V13,V14,V15,V16)
head(result)
result1<-as.data.frame((result))
colnames(result1)<-c("G4_id","name","G4_ref","G4_alt","V9","V10","V11","V12","V13","V14","V15","V16")
##print(head(result1))

###outfile<-paste0(path,"/",out,".RDS")
return(result1)
###saveRDS(g4hunter.input1.2,file=outfile)

}



#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------


myOuterFunction1 <- function(list.x,y) {
    path=dirname("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/coding/coding_indelfiltered.recode.rds")
    tmp.1 <-  dog4hunter(list.x,y,path)
	print("yesssssssssssssssssss")
	return(tmp.1)

}

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

#```
seq.list <- split(codingseq, 1:NROW(codingseq))
codingpath<-dirname("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/coding/coding_indelfiltered.recode.rds")

print("split dataframe")
system.time(result1<-mclapply(seq.list,myOuterFunction1,y="1000genome",mc.cores=8))
###system.time(result1<-lapply(seq.list,myOuterFunction1))

####system.time(result1<-parallel::mclapply(seq.list,dog4hunter(x,y,path), y = "coding",path=codingpath, mc.cores=8))

result1_df<- rbindlist(result1)

####result1_df<- as.data.frame(do.call(rbind, result1))
colnames(result1_df)<-c("G4_id","name","G4_ref","G4_alt","V9","V10","V11","V12","V13","V14","V15","V16")
result1_df$V16<-as.numeric(result1_df$V16)
###result1_df[,5:12] <- sapply(result1_df[,5:12],as.numeric)
print(paste0("TOTAL G4 in CODING:  ", length(which(!is.na(result1_df$V10)))))
print(paste0("VARIANTS WITH NO G4 in CODING:  ", length(which(is.na(result1_df$V10)))))


out<-"coding"

outfile<-paste0(codingpath,"/",out,"_",numrow,".RDS")


saveRDS(result1_df,file=outfile)


print("one completed")
rm(codingseq,result1_df,result1)

end_time <- Sys.time()
difftime(end_time, start_time, units = "secs")[[1]]

exit <- function() { invokeRestart("abort") }    

print("this is the last message")
exit()
print("you should not see this")
#```

ncvpath<-dirname("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/NCV/NCV_indelfiltered.recode.rds")

ncvseq<-readRDS("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/NCV/NCV_indelfiltered.recode.rds")
ncvseq=ncvseq[sample(nrow(ncvseq),replace=F,size=numrow*nrow(ncvseq)),]


print("STARTED NONCODING")

myOuterFunction2 <- function(list.x) {
    y="NCV"
    path=dirname("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/NCV/NCV_indelfiltered.recode.rds")
   
    tmp.1 <-  dog4hunter(list.x,y,path)
	print("yesssssssssssssssssss")
	return(tmp.1)

}

seq.list <- split(ncvseq, 1:NROW(ncvseq))
print("GOING AGAIN")


system.time(result2<-mclapply(seq.list,myOuterFunction2,mc.cores=8 ))
result2_df<- rbindlist(result2)



#result2_df<- data.frame(do.call(rbind, result2))
 colnames(result2_df)<-c("G4_id","name","G4_ref","G4_alt","V9","V10","V11","V12","V13","V14","V15","V16")
result2_df$V16<-as.numeric(result2_df$V16)
#result2_df[,5:12] <- sapply(result2_df[,4:11],as.numeric)
print(paste0("TOTAL G4 in NCV:  ", length(which(!is.na(result2_df$V10)))))
print(paste0("VARIANTS WITH NO G4 in NCV:  ", length(which(is.na(result2_df$V10)))))
out<-"NCV"

end_time <- Sys.time()
difftime(end_time, start_time, units = "secs")[[1]]
outfile<-paste0(ncvpath,"/",out,"_",numrow,".RDS")


saveRDS(result2_df,file=outfile)


print("total done")



myOuterFunction2 <- function(list.x) {
    y="clinvar"
    path=dirname("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/clinvar.rds")
   
    tmp.1 <-  dog4hunter(list.x,y,path)
	print("yesssssssssssssssssss")
	return(tmp.1)
}



clinvarpath<-dirname("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/clinvar.rds")


#clinvarseq=clinvarseq[sample(nrow(clinvarseq),replace=F,size=numrow*nrow(clinvarseq)),]
clinvarseq
seq.list <- split(clinvarseq, 1:NROW(clinvarseq))
print("GOING AGAIN")

system.time(result2<-lapply(seq.list,myOuterFunction2))#,mc.cores=8 ))

system.time(result2<-mclapply(seq.list,myOuterFunction2,mc.cores=8 ))
result2_df<- rbindlist(result2)



#result2_df<- data.frame(do.call(rbind, result2))
 colnames(result2_df)<-c("G4_id","name","G4_ref","G4_alt","V9","V10","V11","V12","V13","V14","V15","V16")
result2_df$V16<-as.numeric(result2_df$V16)
#result2_df[,5:12] <- sapply(result2_df[,4:11],as.numeric)
print(paste0("TOTAL G4 in clinvar:  ", length(which(!is.na(result2_df$V10)))))
print(paste0("VARIANTS WITH NO G4 in clinvar:  ", length(which(is.na(result2_df$V10)))))
out<-"clinvar"

end_time <- Sys.time()
difftime(end_time, start_time, units = "secs")[[1]]
outfile<-paste0(clinvarpath,"/",out,"_",numrow,".RDS")


saveRDS(result2_df,file=outfile)


print("total done")
