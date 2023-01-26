#!/usr/bin/env Rscript
source("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/function_G4Hunter.r")

myOuterFunction1 <- function(list.x,y,path) {
  path=dirname(path_getseq)
  tmp.1 <-  dog4hunter(list.x,y,path)
  i=1
  print(paste0("count = ",i))
  i=i+1	
  return(tmp.1)
  
}
#after COSMIC_getseq,
#input path as required
path_getseq<-"/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/1000gen/vcf1000genome/vcf1000genome.rds"
thoupath<-dirname(path_getseq)

out<-"1000genome"

outfile<-paste0(thoupath,"/",out,".RDS")


thou<-readRDS(path_getseq)

print("read sequence")


seq.list <- split(thou, 1:NROW(thou))

print("split dataframe")
system.time(result1<-mclapply(seq.list,myOuterFunction1,y="1000genome",mc.cores=8))

system.time(result1<-lapply(seq.list,myOuterFunction1,y="1000genome",path=path_getseq))

result1_df<- rbindlist(result1)

####result1_df<- as.data.frame(do.call(rbind, result1))
colnames(result1_df)<-c("G4_id","name","G4_ref","G4_alt","V9","V10","V11","V12","V13","V14","V15","V16")
result1_df$V16<-as.numeric(result1_df$V16)
###result1_df[,5:12] <- sapply(result1_df[,5:12],as.numeric)
print(paste0("TOTAL G4 :  ", length(which(!is.na(result1_df$V10)))))
print(paste0("VARIANTS WITH NO G4 :  ", length(which(is.na(result1_df$V10)))))


out<-"1000genome"


saveRDS(result1_df,file=outfile)



#--------------------------------------------


###result1
count_1000gen<-result1_df[, .(G4hunter_ref = mean(V9), G4hunter_alt = mean(V13),count=.N,sd=sd(V9),sd2=sd(V13)), by = .(V12,V16)]

gen1000_total<-left_join(result1_df,thou,by=c("name","G4_ref","G4_alt","G4_id"))

gen1000_total<-gen1000_total%>% filter(V12==1 | V12==-1 | V16==1 | V16== -1)
gen1000_total<-setDT(gen1000_total)
gen1000_total[, G4_refrc := ifelse((V12 == -1 | V16 == -1) , sapply( G4_ref, function(x) as.character(
  reverseComplement(DNAString(x)))), sapply(G4_ref,function(x) as.character(x)))]

gen1000_total[, G4_altrc:= ifelse(V12 == -1 | V16 == -1, sapply( G4_alt, function(x) as.character(
  reverseComplement(DNAString(x)))), sapply(G4_alt,function(x) as.character(x) ))]
gen1000_total$nameimp<-gen1000_total$name
gen1000_total$name<-gsub("_G4_id","",gen1000_total$name)
gen1000_total<-separate(data = gen1000_total, col = name, into = c("left", "right"), sep = "_")


gen1000_total$start<-as.numeric(gen1000_total$right)-30
gen1000_total$end<-as.numeric(gen1000_total$right)+30

coding2<-gen1000_total
coding2<-coding2%>% filter(nchar(REF)==1 &  nchar(ALT)==1)
coding2$name<-paste0(coding2$left,":",coding2$start,"-",coding2$end,"||","loc:",coding2$loc_start,"||",coding2$REF,"&",coding2$ALT)

#get variants across two strands
table(coding2$V12,coding2$V16)
writeFastaref(coding2,"genome1000_REF.fasta")
writeFastaalt(coding2,"genome1000_ALT.fasta")
#setwd("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/1000gen/")

./../seqquadparser_wrapper.pl /bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/1000gen/genome1000_ALT.fasta \
/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/1000gen/genome1000_ALT_quadparser.tsv

./../seqVariants_cosmic/quadparser_wrapper.pl /bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/1000gen/genome1000_REF.fasta \
/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/1000gen/genome1000_REF_quadparser.tsv

temp_coding_r<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/1000gen/genome1000_REF_quadparser.tsv",sep="\t")
temp_coding_a<-fread("/bio/home/goonerrn/g_quadruplex/annotation_database/KAVIAR/COSMIC/1000gen/genome1000_ALT_quadparser.tsv",sep="\t")
colnames(temp_coding_a)[c(2,3,4,5)]<-c("g4_start","g4_stop","g4seq","g4_type")
colnames(temp_coding_r)[c(2,3,4,5)]<-c("g4_start","g4_stop","g4seq","g4_type")

temp_coding_a%>% group_by(V1)%>% mutate(n=n())%>% arrange(desc(n))%>% head()
temp_coding_a<-temp_coding_a[!duplicated(temp_coding_a), ]
temp_coding_r<-temp_coding_r[!duplicated(temp_coding_r),]

data1<-left_join(coding2,temp_coding_r[,c("V1","g4_start","g4_stop","g4seq","g4_type")],by=c("name"="V1"))
data1na<-data1%>% filter(is.na(g4_start))%>% select(-c(g4_start,g4_stop,g4seq,g4_type))
data1<-data1%>% filter(!is.na(g4_start))


data1na<-left_join(data1na,temp_coding_a[,c("V1","g4_start","g4_stop","g4seq","g4_type")],by=c("name"="V1"))

data_total<-rbind(data1,data1na)



data_total<- data_total%>% 
  mutate(SNP_strand = case_when((V10 == 1 | V16 == 1) ~ merged$SNP,
                                (V10 == -1 | V16 == -1) ~ chartr('ATGC', 'TACG',merged$SNP) ,
                                (V10 == 1 & V16 == -1) ~ merged$SNP,
                                (V10 == -1 & V16 == 1) ~ chartr('ATGC', 'TACG',merged$SNP))) 



data_total<-data_total%>% 
  mutate(G_quadstrand = case_when((V10 == 1 | V16 == 1) ~ "+",
                                  (V10 == -1 | V16 == -1) ~ "-" ,
                                  (V10 == 1 & V16 == -1) ~ "+",
                                  (V10 == -1 & V16 == 1) ~ "-")) 


data_total_1<-data_total%>% as.data.frame()%>%group_by(name,G4_refrc,G4_altrc,V12,V16,start,G_quadstrand,SNP_strand,g4_type)%>% 
  mutate(numberofg4=n(),which1=paste0(unique(which),collapse="|"),g4seq1=paste0(unique(g4seq),collapse="|"),g4_start=min(as.numeric(g4_start)),g4_stop=max(as.numeric(g4_stop)))%>%
  arrange(desc(numberofg4)) %>%
  mutate(g4seq=g4seq1,g4seq1=NULL,which=which1,which=NULL)%>%distinct()%>% as.data.frame()