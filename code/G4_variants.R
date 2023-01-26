library(stringr)

library(data.table)
library(purrr)
library(tidyr)

load("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/g4_G4grindr.RData")
output_fwd<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/query/output_fwd_withCHR.txt",header=FALSE)

g4_exprvalid<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/epimap/additional_bed/pG4_merged_with_GSM3003539_w15.K.merged.sorted.hg19.bed")

a4<-a3%>% mutate(pmaxx=pmax(G4Hunter,G4Hunter_s))%>% group_by(V6)%>% slice(which.max(pmaxx))%>% as.data.frame()

a4$pmaxx_len<-ifelse(a4$Score_s>a4$Score,a4$Length_s,a4$Length)

a4$seq<-ifelse(a4$Score_s>a4$Score,a4$Sequence,a4$V8)

#filter by g4 score
a4<-a4%>% filter(pmaxx>23)

##this is for output_fwd
output_fwd<-output_fwd[output_fwd$V1 %in% chromosome]
output_fwd$name<-paste0(output_fwd$V1,":",output_fwd$V2)
output_fwd<-output_fwd%>%separate(V3,c("nggg","nquad","fquad"),":")
output_fwd$nggg<-as.numeric(output_fwd$nggg)
colnames(output_fwd)[6]<-"sequence"


g4_db<-merge(output_fwd,a4[,c("V6","pmaxx","pmaxx_len","V8")],by.x="name",by.y="V6",all.x=TRUE)

df<-setDT(g4_db)[, pmaxx1 := unique(pmaxx[!is.na(pmaxx)]), by = V4]

df$pmaxx1 <- ave(g4_db$pmaxx, df$V4, FUN=function(x)unique(x[!is.na(x)]))

g4_db.1<-g4_db%>% group_by(V4)%>% mutate(pmaxx1=zoo:na.locf(pmaxx))

g4_db.1%>% is.na(pmaxx1)%>%nrow()
g4_db%>% group_by(V4) %>%  filter(all(is.na(pmaxx)))%>% nrow()



g4_expr<-merge(db_db,g4_exprvalid,by.x="name",by.y="V5")'
here, the V1.y, V2.y, V3 are hg19 positions.

g4_expr%>% is.na(pmaxx1)%>%nrow()

#we make hg38 bed file to look for variants in COSMIC VCF files,
g4_expr<-separate(g4_expr,name,into=c("chr","range"),sep=":",remove=FALSE)
g4_expr<-separate(g4_expr,range,into=c("start","end"),sep="-",remove=TRUE)
g4_expr.1<-g4_expr%>% select(chr,start,end, name, pmaxx)
g4_expr.1$pmaxx[is.na(g4_expr.1$pmaxx)]

g4_expr.1$chr<-gsub("chr","",g4_expr.1$chr)

fwrite(g4_expr.1, "/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/g4_expr.bed",sep='\t',col.names=FALSE,quote=FALSE)







#----------------------------------------------------------
##get the bed file intersected with variants and then
##either or

library(dplyr)
g4_var<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/g4_exp_NCV.bed")%>% distinct()
g4_var<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/g4_codingvar.bed")%>% distinct()


row_tibble <- function(x, col_names) {
  tibble::as_tibble(rbind(setNames(x, col_names)))
}

parse_info <- function(info) {
  strsplit(info, ";", fixed = TRUE) %>%
    purrr::map(~row_tibble(sub("^.*=(.*)", "\\1", .x), sub("^(.*)=.*", "\\1", .x)))%>% bind_rows(.id="source")
}

g4_info_parsed<-parse_info(g4_var$V13)
g4_info_parsed<-separate(g4_info_parsed,GENE,into=c("gene","transcript"),sep="_",remove=TRUE)
g4_var<-cbind(g4_var,g4_info_parsed)

g4_var$V13<-NULL

g4_var<-g4_var%>% group_by(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,STRAND,LEGACY_ID,CNT)%>% unique()%>%
summarise(transcripts=paste0(unique(na.omit(transcript)),collapse=";"),genes=paste0(unique(na.omit(gene)),collapse=";"))%>% as.data.frame()


g4_var$var<-paste0(g4_var$V9,"->",g4_var$V10)

##need the sequence for comparison later
g4_var<-merge(g4_var,output_fwd[,c("name","sequence")],by.x="V4",by.y="name")

#fwrite(g4_var, "/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/g4_var_NCV.tsv",sep='\t',col.names=FALSE,quote=FALSE)
fwrite(g4_var, "/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/g4_var_codingV.tsv",sep='\t',col.names=FALSE,quote=FALSE)


#_-----------------------------------------------------
##annovar

#for coding
library(data.table)
g4_var<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/g4_var_codingV.tsv")%>% distinct()
hg38_coding<-fread("bio/home/goonerrn/tools/annovar/hg38_cosmic95")
hg38_coding$V2<-as.numeric(hg38_coding$V2)
##remove gene to keep one G quad per code, 
g4_var.1<-g4_var%>% select(-V11)%>% distinct()
g4_var_annovar<-inner_join(g4_var.1,hg38_coding,by=c("V6"="V1","V7"="V2"))


g4_info_parsed<-parse_info(g4_var_annovar$V6.y)
g4_info_parsed<-g4_info_parsed %>% 
    mutate(OCCURENCE = strsplit(as.character(OCCURENCE), ",")) %>% 
    unnest(OCCURENCE)
g4_info_parsed<-separate(g4_info_parsed,OCCURENCE,into=c("NUMBER","CELL_TYPE"),sep="[()]",remove=TRUE)


g4_var.2<-inner_join(g4_var.1,g4_info_parsed,by=c("V8"="ID"))

#g4_var$V13<-NULL

#----------------------

#for non coding

mutantexport<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/MutantExport.tsv")


g4_var<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/g4_codingNCV.tsv")


left_join(g4_var,mutantexport,by=("GENOMIC_MUTATION_ID"="V2"))
temp<-left_join(g4_var,mutantexport,by=c("V8"="GENOMIC_MUTATION_ID"))
temp<-merge(g4_var,mutantexport,by.x="V8",by.y="GENOMIC_MUTATION_ID")
