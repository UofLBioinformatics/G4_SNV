library(tibble)
library(dplyr)
library(tidyr)

library(dplyr)
library(purrr)
library(ggplot2)

require(Biostrings) 
library(data.table)

writeFastaref<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"G4_refrc"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}


writeFastaalt<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"G4_altrc"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

coding<-readRDS("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/coding/coding_1.RDS")
clinvar_0<-readRDS("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/clinvar.RDS")

 clinvar_2<-left_join(clinvar,clinvar_0[,c("name","REF","ALT")])

coding1<-coding
rm(coding)
coding1<-coding1%>% filter(V12==1 | V12==-1 | V16==1 | V16== -1)
coding1<-setDT(coding1)
coding1[, G4_refrc := ifelse((V12 == -1 | V16 == -1) , sapply( G4_ref, function(x) as.character(
reverseComplement(DNAString(x)))), sapply(G4_ref,function(x) as.character(x)))]


coding1[, G4_altrc:= ifelse(V12 == -1 | V16 == -1, sapply( G4_alt, function(x) as.character(
reverseComplement(DNAString(x)))), sapply(G4_alt,function(x) as.character(x) ))]

coding1$name<-gsub("*_G4_id","",coding1$name)
coding1<-separate(data = coding1, col = name, into = c("left", "right"), sep = "_")

coding1$start<-as.numeric(coding1$right)-30
coding1$end<-as.numeric(coding1$right)+30
coding2$name<-paste0(coding2$left,":",coding2$start,"-",coding2$end,"||","loc:",coding2$loc_start,"||",coding2$REF,"&",coding2$ALT)
coding2<-coding1
writeFastaref(coding2,"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/coding/coding_1_ref_g4.fasta")
writeFastaalt(coding2,"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/coding/coding_1_alt_g4.fasta")
#this is repeat for clinvar
#writeFastaref(coding2,"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/clinvar_1_ref_g4.fasta")
#writeFastaalt(coding2,"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/clinvar/clinvar/clinvar_1_alt_g4.fasta")

print(nrow(coding2))
print(nrow(coding1))
rm(coding1)





myOuterFunction2 <- function(list.x) {
    y="NCV"
    path=dirname("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/seqVariants_cosmic/NCV/NCV_indelfiltered.recode.rds")
    tmp.1 <-  dog4hunter(list.x,y,path)
	print("yesssssssssssssssssss")
	return(tmp.1)

}

