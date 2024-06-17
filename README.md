#Select the data

library(tidyverse)
SumData <- my_data %>% select
exon2 <- exon %>% select ("##date: 2023-09-19","V4","V5")
EUR2 <- EUR1 %>% select(4,5,6)



#creat GRange file 

makeGRangesFromDataFrame(data.frame(SumData),  keep.extra.columns=FALSE,  ignore.strand=FALSE, seqinfo=NULL,  seqnames.field=c("seqnames", "seqname", "chromosome", "chrom",  "chr", "chromosome_name", "seqid"), start.field="start", end.field=c("end", "stop"), starts.in.df.are.0based=FALSE, na.rm=TRUE)
pophuman.region <- makeGRangesFromDataFrame(data.frame(SumData),  keep.extra.columns=FALSE,  ignore.strand=FALSE, seqinfo=NULL,  seqnames.field=c("seqnames", "seqname", "chromosome", "chrom",  "chr", "chromosome_name", "seqid"), start.field="start", end.field=c("end", "stop"), starts.in.df.are.0based=FALSE, na.rm=TRUE)

#change the name of the clome

colnames(RNA.anotation.summery) <- c("chr","start","end")

#Overlab 

numOverlaps(pophuman.region, RNA.anotationGrange, count.once=TRUE)
pt <- overlapPermTest(A=pophuman.region, B=RNA.anotationGrange, ntimes=50)
numOverlaps(randomizeRegions(pophuman.region), RNA.anotationGrange, count.once=TRUE)
numOverlaps(AFRGrange, longnoncodingGrange, count.once=TRUE)
ptAFR <- overlapPermTest(A=AFRGrange, B=longnoncodingGrange, ntimes=100)
lz <- localZScore(pt=pt, A=pophuman.region, B=transcript3)
plot(lz)


#select data acordin row name

gene <- RNA_anotation[RNA_anotation$V3 %in% c('gene'), ]
  
#filter rows based on colounm

EUR <- na.omit(listRegions[, "EUR"])
EUR1 <- listRegions[!(is.na(listRegions$EUR) | listRegions$EUR == ""), ]
filtered_data <- listRegions[listRegions$SFS != 0, ]


#All

Recent1 <- listRegions[!(is.na(listRegions$Recent) | listRegions$Recent == ""), ]
Recent2 <- Recent1 %>% select(4,5,6)
RecentGrang <- makeGRangesFromDataFrame(data.frame(Recent2),  keep.extra.columns=FALSE,  ignore.strand=FALSE, seqinfo=NULL,  seqnames.field=c("seqnames", "seqname", "chromosome", "chrom",  "chr", "chromosome_name", "seqid"), start.field="start", end.field=c("end", "stop"), starts.in.df.are.0based=FALSE, na.rm=TRUE)
numOverlaps(RecentGrang, longnoncodingGrange, count.once=TRUE)
ptRecent <- overlapPermTest(A=RecentGrang, B=longnoncodingGrange, ntimes=100)
ptRecent
plot(ptRecent)
