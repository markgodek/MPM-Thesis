library(ggplot2)
library(ggrepel)
library(data.table)
library(cowplot)
library(stringr)
theme_set(theme_cowplot())

AD_cutoff <- 15 #sets cutoff for insufficient allele depth
DP_Cutoff <- 15 #sets cutoff for insufficient read depth 
plot_height<-5
plot_width<-6
geom_size<-3

setwd('C:/Users/Mark/Dropbox (Partners HealthCare)/Mark Shared MPM Folder') #home 
#setwd('C:/Users/mi313/Dropbox (Partners HealthCare)/Mark Shared MPM Folder') #work

process.cosmic.tsv <- function(cosmic.file){
  MutantColumns<-fread(file=cosmic.file)
  
  MutantColumns<-MutantColumns[SNP=='y',]                                                 #get rows where the mutation is a SNP
  names(MutantColumns)[names(MutantColumns)=="Mutation genome position"] <- "interval"    #rename a column
  
  #pull out CHROM, POS, REF, and ALT from the interval and remove columns we don't need
  MutantColumns[,CHROM:=gsub(":.*","",interval)]
  MutantColumns[,POS:=gsub(".*-","",interval)]
  MutantColumns[,REF:=gsub('c.*[0-9]',"", `Mutation CDS`)]
  MutantColumns[,REF:=gsub('>.',"", REF)]
  MutantColumns[,ALT:=gsub('.*>',"", `Mutation CDS`)]
  MutantColumns<-MutantColumns[,-c('interval','Mutation CDS')]
  
  #divide the dataset into rows which had duplicates and unique rows, i.e. dup.MutantColumns is a subset of MutantColumns - 681051 vs 727231
  dup.MutantColumns<-distinct(dup.MutantColumns[order(CHROM,POS,ALT)])
  MutantColumns<- distinct(MutantColumns,CHROM,POS,ALT, .keep_all=TRUE)
  
  write.table(MutantColumns, file='COSMIC-snp.no-dups.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
  write.table(dup.MutantColumns, file='COSMIC-snp.dups.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
  
}

#process.cosmic.tsv('../CosmicColumns.w.ALT.tsv')

#read the germline alleles and melt the wide format into long format based on the columns listed
af.dt<-fread(file="AF/germline/germline.PASS.AFonly.frq") 
af.melt.dt<-melt(af.dt, id.vars = c("CHROM", "POS", "REF", "ALT", 
                                    "AF", "GNOMAD_AF",  "GENE",  "FEATURE", "DP", "AC", "AN"))

#create a dataframe for the genotype data and extract which place the sample name in a new column
af.GT.melt.dt<-af.melt.dt[variable%like%"GT",]
colnames(af.GT.melt.dt)[13]<-"genotype"
af.GT.melt.dt[,sample:=gsub(":.*","",variable)]

#create a dataframe for the allele depth data, extract the REF counts, ALT counts, and sample name and place them in new columns
af.AD.melt.dt<-af.melt.dt[variable%like%"AD",]
af.AD.melt.dt[,refCount:=as.numeric(gsub(",.*","",value)) ]
af.AD.melt.dt[,altCount:=as.numeric(gsub(".*,","",value)) ]
af.AD.melt.dt[,sample:=gsub(":.*","",variable)]

keycols <- c('CHROM', 'POS', 'REF', 'ALT', 'AF', 
             'GNOMAD_AF', 'GENE', 'FEATURE', 'DP', 'AC', 'AN', 'sample')

#merge the allele depth and genotype dataframes on the key columns and melt them into the final dataframe
setkeyv(af.GT.melt.dt, keycols)
setkeyv(af.AD.melt.dt, keycols)
final.frame <- af.GT.melt.dt[af.AD.melt.dt, ]

remove(af.dt,af.melt.dt,af.AD.melt.dt,af.GT.melt.dt)

colnames(final.frame)[16]<-'Counts'
final.frame<-final.frame[,.(CHROM, POS, REF, ALT, AF, GNOMAD_AF, GENE,
                            FEATURE, DP, AC, AN, sample, genotype,Counts, refCount, altCount)]

final.frame[,genotype_final:=gsub("\\|","/",genotype)]                            #use the genotype to make a formatted genotype column called "genotype_final"
final.frame[refCount<AD_cutoff, genotype_final:=gsub("0","i",genotype_final) ]    #if the REF or ALT allele counts are below a cutoff value, replace them with "i"
final.frame[altCount<AD_cutoff, genotype_final:=gsub("1","i",genotype_final) ]
final.frame[,genotype_final:=gsub("i/1","1/i",genotype_final) ]                   #convert "1/i" genotypes to "i/1"

#store the original genotype, invalidate genotypes with counts below our threshold by assigning them "i" values, and overwrite the genotype_final column
final.frame[,genotype_old:=genotype_final]
final.frame[refCount<AD_cutoff&altCount<AD_cutoff,genotype_final:='i/i']
final.frame[refCount<AD_cutoff&altCount>=AD_cutoff,genotype_final:='i/1']
final.frame[refCount>=AD_cutoff&altCount<AD_cutoff,genotype_final:='0/i']
final.frame[refCount>=AD_cutoff&altCount==0,genotype_final:='0/0']
final.frame[refCount==0&altCount>=AD_cutoff,genotype_final:='1/1']
final.frame[refCount>=AD_cutoff&altCount>=AD_cutoff,genotype_final:='0/1']
final.frame[genotype=='./.', genotype_final:='./.']

#if the row is a member of one of the families, store the sample name in a new column
family.A <- c('HITS622847','HITS622849')
family.B <- c('HITS622856','HITS622858','HITS622870')

#for each CHROM and POS, extract genotypes if the samples are a member of a family and store it in a new column
final.frame[sample %in% family.A,family.A.GT:=paste(genotype_final, sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT")]
final.frame[sample %in% family.B,family.B.GT:=paste(genotype_final, sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT")]

#create new dataframes with only members of each family and create a new
A.frame<-unique(final.frame[family.A.GT!='<NA>',c('CHROM','POS','REF','ALT','family.A.GT')])
B.frame<-unique(final.frame[family.B.GT!='<NA>',c('CHROM','POS','REF','ALT','family.B.GT')])
keycols<-c("CHROM","POS","REF","ALT")
setkeyv(A.frame,keycols)
setkeyv(B.frame,keycols)
family.frame<-merge(x = A.frame, y = B.frame, by = keycols, all.x = TRUE)

#store old AC and AN, then recalculate based off genotype_final which is determined for our cutoffs for reassigning "i" genotypes
final.frame[,AC_old:=AC]
final.frame[,AN_old:=AN]
final.frame[,AC_COMP:=ifelse(genotype_final=='1/1',2,ifelse(genotype_final%in%c('i/1','0/1'),1,0))]
final.frame[,AN_COMP:=ifelse(genotype_final%in%c('1/1','0/0','0/1'),2,ifelse(genotype_final%in%c('i/1','0/i'),1,0))]
final.frame[,AC:=sum(AC_COMP),by=keycols]
final.frame[,AN:=sum(AN_COMP),by=keycols]

#choose columns for final data table and collapse samples from long format back to wide format
final.dt<-final.frame[,.N,by=c('CHROM', 'POS', 'REF', 'ALT', 'AF', 'GNOMAD_AF', 'GENE', 'FEATURE', 'DP', 'AC', 'AN','AC_old','AN_old','genotype_final')]
final.dt<-dcast(final.dt, CHROM+POS+REF+ALT+AF+GNOMAD_AF+GENE+FEATURE+DP+AC+AN+AC_old+AN_old~genotype_final,value.var='N', fill=0)

#store old allele frequency and calculate a new one based on the updated AC and AN
final.dt[,AF_old:=AF]
final.dt[,AF:=AC/AN]
final.dt<-final.dt[AC>0&AN>0,]

#add the family genotype data to the final data table
setkeyv(family.frame,keycols)
setkeyv(final.dt,keycols)
final.dt<-merge(x=final.dt,y=family.frame,by=keycols,all.x = TRUE)

#read COSMIC data and add it to our table based on CHROM, POS, and ALT
cosmic.dt<-fread(file = 'COSMIC-snp.no-dups.tsv')
cosmic.dt$CHROM<-as.character(cosmic.dt$CHROM)
keycols<-c("CHROM","POS","ALT")
setkeyv(final.dt,keycols)
setkeyv(cosmic.dt,keycols)
final.w.cosmic.ALT<-merge(x = final.dt, y = cosmic.dt, by = keycols, all.x = TRUE)

remove(final.frame,final.dt,A.frame,B.frame,family.frame,cosmic.dt)

#read the CADD data, rename some columns, and add it to our table based on CHROM, POS, and ALT
CADD.dt<-fread(file = 'CADD_annotations_for_germline.PASS.tsv')
colnames(CADD.dt)[1]<-'CHROM'
colnames(CADD.dt)[2]<-'POS'
colnames(CADD.dt)[3]<-'REF'
colnames(CADD.dt)[4]<-'ALT'
setkeyv(final.w.cosmic.ALT,keycols)
setkeyv(CADD.dt,keycols)
final.w.CADD.cosmic<-merge(x = final.w.cosmic.ALT, y = CADD.dt, by = keycols, all.x = TRUE)
final.w.CADD.cosmic<-final.w.CADD.cosmic[,!c("REF.y","REF.x")]

#combine the columns below for each set of unique values of CHROM, POS, REF, and ALT so the rows will collapse nicely without losing data
final.w.CADD.cosmic[,Consequence:=paste(unique(Consequence), sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT")]
final.w.CADD.cosmic[,AnnoType:=paste(AnnoType, sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT")]
final.w.CADD.cosmic[,motifEName:=paste(motifEName, sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT")]
final.w.CADD.cosmic[,GeneID:=paste(GeneID, sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT")]
final.w.CADD.cosmic[,FeatureID:=paste(FeatureID, sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT")]
final.w.CADD.cosmic[,GeneName:=paste(GeneName, sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT")]
final.w.CADD.cosmic[,MaxConsScore:=paste(max(ConsScore, sep = "", collapse = ", ")),by=c("CHROM","POS","REF","ALT")]
final.w.CADD.cosmic[MaxConsScore=='NA', MaxConsScore:=0]

#select columns and remove duplicate rows
final.w.CADD.cosmic<-unique(final.w.CADD.cosmic[,c('CHROM','POS', 'REF', 'ALT','AF','AF_old','GNOMAD_AF', 'GENE','FEATURE','DP',
                                            'AC','AN','AC_old','AN_old','./.','0/0','0/1','0/i','1/1','i/1','i/i','family.A.GT',
                                            'family.B.GT','FATHMM prediction', 'Mutation somatic status', 'Consequence', 'GeneID','MaxConsScore','PHRED')])

remove(final.w.cosmic.ALT,CADD.dt)

write.table(final.w.CADD.cosmic, file="final.variants.annotated.tsv", sep="\t", row.names=FALSE)

#add columns for genotypes of families
final.w.CADD.cosmic[,"A.0/0":=str_count(family.A.GT,'0/0')]
final.w.CADD.cosmic[,"A.0/i":=str_count(family.A.GT,'0/i')]
final.w.CADD.cosmic[,"A.0/1":=str_count(family.A.GT,'0/1')]
final.w.CADD.cosmic[,"A.1/1":=str_count(family.A.GT,'1/1')]
final.w.CADD.cosmic[,"A.i/1":=str_count(family.A.GT,'i/1')]
final.w.CADD.cosmic[,"A.i/i":=str_count(family.A.GT,'i/i')]
final.w.CADD.cosmic[,"A../.":=str_count(family.A.GT,'\\./\\.')]
final.w.CADD.cosmic[,"B.0/0":=str_count(family.B.GT,'0/0')]
final.w.CADD.cosmic[,"B.0/1":=str_count(family.B.GT,'0/1')]
final.w.CADD.cosmic[,"B.0/i":=str_count(family.B.GT,'0/i')]
final.w.CADD.cosmic[,"B.1/1":=str_count(family.B.GT,'1/1')]
final.w.CADD.cosmic[,"B.i/1":=str_count(family.B.GT,'i/1')]
final.w.CADD.cosmic[,"B.i/i":=str_count(family.B.GT,'i/i')]
final.w.CADD.cosmic[,"B../.":=str_count(family.B.GT,'\\./\\.')]

#count number of zero or 1. if there are 2 or more i's, reduce observations to 1, else either two 1's or 0/1 were seen, i.e. 2 observations
#add these observations from AN
final.w.CADD.cosmic[,'AN+.A':=ifelse(`A.1/1`>=1|`A.0/1`>=1|`A.0/0`>=1,2,ifelse(`A.i/1`>=1|`A.0/1`>=1,1,0))]
final.w.CADD.cosmic[,'AN+.B':=ifelse(`B.1/1`>=1|`B.0/1`>=1|`B.0/0`>=1,2,ifelse(`B.i/1`>=1|`B.0/1`>=1,1,0))]

#subtract zero and 1 observations from AN and 1's from AC
#if there is a genotype, count zeroes and 1's, else set it to zero - code to deal with NA's in genotype column
final.w.CADD.cosmic[,'0.counts':=str_count(family.A.GT,'0')+str_count(family.B.GT,'0')]
final.w.CADD.cosmic[,'1.counts':=str_count(family.A.GT,'1')+str_count(family.B.GT,'1')]

#count number of 1's in the genotypes. 1/1 = 2 and overides 0/1 and i/1, which = 1
final.w.CADD.cosmic[,'AC+.A':=ifelse(`A.1/1`>=1,2,ifelse(`A.0/1`>=1|`A.i/1`>=1,1,0))]
final.w.CADD.cosmic[,'AC+.B':=ifelse(`B.1/1`>=1,2,ifelse(`B.0/1`>=1|`B.i/1`>=1,1,0))]

#adjusted AN - AN with observations in all families flattened
final.w.CADD.cosmic[,'adj_AN':=AN+`AN+.A`+`AN+.B`-`1.counts`-`0.counts`]

#adjusted AC - AC with all observations in families flattened
final.w.CADD.cosmic[,'adj_AC':=AC-`1.counts`+`AC+.A`+`AC+.B`]

#adjusted AF - AF of variant after reducing observations per family to one
final.w.CADD.cosmic[,adj_AF:= adj_AC/adj_AN]

#add some columns useful for filtering - high AF_distance means 
final.w.CADD.cosmic[,'family':=`A.1/1`+`A.0/1`+`A.i/1`+`B.1/1`+`B.0/1`+`B.i/1`]
final.w.CADD.cosmic[,'AF_dist':=abs(adj_AF-as.numeric(GNOMAD_AF))]
final.w.CADD.cosmic[,'LRFR':=log10(adj_AF/as.numeric(GNOMAD_AF))]
final.w.CADD.cosmic[,'Rare':=ifelse(LRFR>0.75|GNOMAD_AF=='.','Yes','No')]       # a variant is "rare" if it has a high log ratio, or if it does not appear in the gnomAD dataset
final.w.CADD.cosmic[,'Impactful':=ifelse(is.na(PHRED),'No',ifelse(PHRED>=20,'Yes','No'))]
final.w.CADD.cosmic[,'Important':=ifelse((Rare=='Yes'&Impactful=='Yes'),'Yes','No')]
final.w.CADD.cosmic[(GNOMAD_AF=='.'&PHRED>=20&adj_AF>0.125), 'Important':='Yes']         # Two samples that are important do not have GNOMAD_AF values, so annotated them manually

#create a column to color the points in the plots with the heirarch important > impactful > rare
#final.w.CADD.cosmic[,Color_explicit:=ifelse(Important=='Yes','red',ifelse(Impactful=='Yes','orange',ifelse(Rare=='Yes','gold','black')))]
final.w.CADD.cosmic[,Color:=ifelse(Important=='Yes','Important',ifelse(Impactful=='Yes','Impactful',ifelse(Rare=='Yes','Rare','Common')))]

#extract different permutations of family members into new data tables
familyA.dt<-final.w.CADD.cosmic[(`A.1/1`+`A.0/1`+`A.i/1`)==2]                 #familyA with 2 members
familyB2.dt<-final.w.CADD.cosmic[(`B.1/1`+`B.0/1`+`B.i/1`)==2]                                             #familyB with 2 members
familyB3.dt<-final.w.CADD.cosmic[(`B.1/1`+`B.0/1`+`B.i/1`)==3]                #familyB with 3 members
all.familyB.dt<-rbind(familyB2.dt,familyB3.dt)                                                  #variants present in 2 or 3 members of family B
any.siblings.dt<-rbind(all.familyB.dt,familyA.dt)                                               #variants present in both families with 2 or more members
all.five.dt<-final.w.CADD.cosmic[(`A.1/1`+`A.0/1`+`A.i/1`+`B.1/1`+`B.0/1`+`B.i/1`)==5]                     #variants present in all 5 members of the families
only.family.all.five.dt<-final.w.CADD.cosmic[(`0/1`+`1/1`+`i/1`)<=10&(`A.1/1`+`A.0/1`+`A.i/1`+`B.1/1`+`B.0/1`+`B.i/1`)==5]

#extract data tables with were there are only variants within family A or B
only.family.A.dt<-fsetdiff(familyA.dt,all.five.dt)
only.family.B.dt<-fsetdiff(familyB3.dt,all.five.dt)

final.no.GNOMAD <-final.w.CADD.cosmic[GNOMAD_AF=='.'] 
final.w.GNOMAD <-final.w.CADD.cosmic[GNOMAD_AF!='.']

keeperColumns<-c('CHROM','POS', 'REF', 'ALT','AF','AF_old','GNOMAD_AF', 'GENE','FEATURE','DP',
                 'AC','AN','AC_old','AN_old','./.','0/0','0/1','0/i','1/1','i/1','i/i','family.A.GT',
                 'family.B.GT','FATHMM prediction', 'Mutation somatic status', 'Consequence', 'GeneID',
                 'MaxConsScore','adj_AN','adj_AC','adj_AF','LRFR', 'PHRED', 'Rare', 'Impactful', 'Important', 'Color')



final.w.CADD.cosmic<-final.w.CADD.cosmic[, ..keeperColumns]
pathogenic.cosmic.dt<-final.w.CADD.cosmic[`FATHMM prediction`=="PATHOGENIC"] #final dataframe for variants pathogenic according to COSMIC data
final.no.GNOMAD<-final.no.GNOMAD[, ..keeperColumns]
final.w.GNOMAD<-final.w.GNOMAD[, ..keeperColumns]
all.five.dt<-all.five.dt[, ..keeperColumns]
only.family.A.dt<-only.family.A.dt[, ..keeperColumns]
only.family.B.dt<-only.family.B.dt[, ..keeperColumns]

no.PHRED.but.Rare<-final.w.CADD.cosmic[is.na(PHRED)][Rare=='Yes']   #report these

# create output directory for plots
plot.path<-'/Results/germline/plots'
dir.create(paste0(getwd(),plot.path), showWarnings = FALSE)

# figure 1.1.1a - gnomad vs cohort frequency to describe what the ratio is - a metric for rarity
ggplot(final.w.CADD.cosmic[DP>DP_Cutoff,][(`0/1`+`1/1`+`i/1`)>0, ],aes(x=adj_AF, y=-log10(as.numeric(GNOMAD_AF))))+
    geom_point(alpha = 0.08,size = geom_size)+geom_hline(yintercept = -log10(0.01), color="red")+
    geom_vline(xintercept = 0.0625, color="red") + 
    labs(x='Adjusted Cohort Allele Frequency',y='gnomAD v2.1 Allele Frequency (-log scale)') +
    theme(plot.title = element_text(hjust = 0.5))
ggsave('Figure 1.1.1a - gnomAD vs AF.png', plot = last_plot(), device = 'png', path = paste0(getwd(),plot.path), scale = 1, height = plot_height, width = plot_width, dpi = 'screen', units = 'in', limitsize = TRUE)    

#only variants which had a gnomAD AF
insetPlot<-final.w.CADD.cosmic[DP>DP_Cutoff,][(`0/1`+`1/1`+`i/1`)>0, ][as.numeric(GNOMAD_AF)<0.01&adj_AF>=0.0625,]
insetPlot[GENE=='Unknown'|GENE%like%'_'|GENE%like%'\\.',GENE:=NA]

# figure 1.1.1b - inset of figure 1a on interesting quadrant
ggplot(insetPlot,aes(x=adj_AF, y=-log10(as.numeric(GNOMAD_AF))))+
    labs(x='Adjusted Cohort Allele Frequency',y='gnomAD v2.1 Allele Frequency (-log scale)') +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_text_repel(data=subset(insetPlot, adj_AF > 0.25), 
    aes(label=GENE)) +
    geom_point(alpha = 0.08,size = geom_size)
ggsave('Figure 1.1.1b - gnomAD vs AF - inset.png', plot = last_plot(), device = 'png', path = paste0(getwd(),plot.path), scale = 1, height = plot_height, width = plot_width, dpi = 'screen', units = 'in', limitsize = TRUE)    

# figure 1.1.2 - log vs absolute difference - showing that we considered different metrics to determine rarity
ggplot(final.w.CADD.cosmic,aes(x=(adj_AF-as.numeric(GNOMAD_AF)), y=LRFR))+
    labs(x='Adjusted Cohort Allele Frequency - gnomAD Allele Frequency',y='Log Relative Frequency Ratio') +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_point(alpha = 0.08,size = geom_size)
ggsave('Figure 1.1.2 - log vs abs diff.png', plot = last_plot(), device = 'png', path = paste0(getwd(),plot.path), scale = 1, height = plot_height, width = plot_width, dpi = 'screen', units = 'in', limitsize = TRUE)    

# figure 1.1.3 - histogram of log relative allele frequency for whole dataset - 0.75 threshold for rarity
ggplot(final.w.CADD.cosmic[GNOMAD_AF!='.',]) +
    aes(LRFR) +
    geom_histogram(binwidth = 0.25,color='black',fill='white') +
    labs(x='Log Relative Frequency Ratio', y='# of Variants') +
    geom_vline(xintercept = 0.75, color="red") + 
    theme(plot.title = element_text(hjust = 0.5), plot.caption = element_text(size = 10))
ggsave('Figure 1.1.3 - AF - hist.png', plot = last_plot(), device = 'png', path = paste0(getwd(),plot.path), scale = 1, height = plot_height, width = plot_width, dpi = 'screen', units = 'in', limitsize = TRUE)    

# figure 1.1.4 - scatter plot of PHRED and Consequence - shows we considered different metrics to determine impact
ggplot(final.w.CADD.cosmic,aes(x=MaxConsScore, y=PHRED))+
    labs(x='Maximum Consequence Score',y='PHRED',title=paste0('')) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_point(alpha = 0.08,size = geom_size)
ggsave('Figure 1.1.4 - MaxCons vs PHRED.png', plot = last_plot(), device = 'png', path = paste0(getwd(),plot.path), scale = 1, height = plot_height, width = plot_width, dpi = 'screen', units = 'in', limitsize = TRUE)    

# figure 1.1.5 - histogram of PHRED score - to determine cutoff for variants of interest - score of greater or equal 20 indicates the 1% most deleterious (https://cadd.gs.washington.edu/info)
ggplot(final.w.CADD.cosmic) +
    aes(PHRED) +
    geom_histogram(binwidth = 1,color='black',fill='white') +
    labs(x='PHRED', y='# of Variants') +
    geom_vline(xintercept = 20, color="red") +
    theme(plot.title = element_text(hjust = 0.5), plot.caption = element_text(size = 10))
ggsave('Figure 1.1.5 - PHRED - hist.png', plot = last_plot(), device = 'png', path = paste0(getwd(),plot.path), scale = 1, height = plot_height, width = plot_width, dpi = 'screen', units = 'in', limitsize = TRUE)    

# figure 1.2.1 - PHRED vs ratio - with GNOMAD - dots red if important - gold if rare - orange if impactful
ggplot(final.w.GNOMAD,aes(x=LRFR, y=PHRED, color=Color))+
    labs(x='Log Relative Frequency Ratio',y='PHRED') +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_point(alpha = 0.08,size = geom_size) + 
    scale_color_manual(name = "",breaks = c('Rare', 'Important','Impactful'), 
                       values=c('black','orange','red','gold'), guide=guide_legend(reverse=TRUE)) +
    guides(color = guide_legend(override.aes = list(alpha = 3) ) )
ggsave('Figure 1.2.1 - Ratio vs PHRED - w gnomAD.png', plot = last_plot(), device = 'png', path = paste0(getwd(),plot.path), scale = 1, height = plot_height, width = plot_width, dpi = 'screen', units = 'in', limitsize = TRUE)    

#figure 1.2.2 - PHRED vs adj_AF - no gnomad - dots red if important - gold if rare - orange if impactful
ggplot(final.no.GNOMAD,aes(x=adj_AF, y=PHRED, color=Color))+
    labs(x='Adjusted Cohort Allele Frequency',y='PHRED') +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_point(alpha = 0.08,size = geom_size) + 
    scale_color_manual(name = "",breaks = c('Rare', 'Important','Impactful'), 
                       values=c('gold','red'), guide=guide_legend(reverse=TRUE)) +
    guides(color = guide_legend(override.aes = list(alpha = 3) ) )
ggsave('Figure 1.2.2 - Adj_AF vs PHRED - no gnomAD.png', plot = last_plot(), device = 'png', path = paste0(getwd(),plot.path), scale = 1, height = plot_height, width = plot_width, dpi = 'screen', units = 'in', limitsize = TRUE)    

#figure 1.3.1a
ggplot(all.five.dt[GNOMAD_AF!='.'],aes(x=LRFR, y=PHRED, color=Color)) +
    labs(x='Log Relative Frequency Ratio',y='PHRED') +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_point(alpha = 0.15,size = geom_size) + 
    scale_color_manual(name = "", breaks = c('Rare', 'Impactful'),
                       values=c('black','orange','gold'), guide=guide_legend(reverse=TRUE)) +
    guides(color = guide_legend(override.aes = list(alpha = 3) ) )
ggsave('Figure 1.3.1a - all five - RLFR vs PHRED.png', plot = last_plot(), device = 'png', path = paste0(getwd(),plot.path), scale = 1, height = plot_height, width = plot_width, dpi = 'screen', units = 'in', limitsize = TRUE)    

#figure 1.3.1b
ggplot(all.five.dt[GNOMAD_AF=='.'],aes(x=adj_AF, y=PHRED, color=Color)) +
    labs(x='Adjusted Cohort Allele Frequency',y='PHRED') +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_point(alpha = 1,size = geom_size) + 
    scale_color_manual(name = "", breaks = c('Rare', 'Important','Impactful'), 
                       values=c('gold'), guide=guide_legend(reverse=TRUE)) +
    guides(color = guide_legend(override.aes = list(alpha = 3) ) )
ggsave('Figure 1.3.1b - all five - Adj_AF vs PHRED.png', plot = last_plot(), device = 'png', path = paste0(getwd(),plot.path), scale = 1, height = plot_height, width = plot_width, dpi = 'screen', units = 'in', limitsize = TRUE)    

#figure 1.3.2a
ggplot(only.family.A.dt[GNOMAD_AF!='.'],aes(x=LRFR, y=PHRED, color=Color)) +
    labs(x='Log Relative Frequency Ratio',y='PHRED') +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_point(alpha = 0.15,size = geom_size) + 
    scale_color_manual(name = "", breaks = c('Rare', 'Important','Impactful'), 
                       values=c('black','orange','red','gold'), guide=guide_legend(reverse=TRUE)) +
    guides(color = guide_legend(override.aes = list(alpha = 3) ) )
ggsave('Figure 1.3.2a - family A - LRFR vs PHRED.png', plot = last_plot(), device = 'png', path = paste0(getwd(),plot.path), scale = 1, height = plot_height, width = plot_width, dpi = 'screen', units = 'in', limitsize = TRUE)    

#figure 1.3.2b
ggplot(only.family.A.dt[GNOMAD_AF=='.'],aes(x=adj_AF, y=PHRED, color=Color)) +
    labs(x='Adjusted Cohort Allele Frequency',y='PHRED') +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_point(alpha = 0.1,size = geom_size) + 
    scale_color_manual(name = "", breaks = c('Rare', 'Important','Impactful'), 
                       values=c('gold','red'), guide=guide_legend(reverse=TRUE)) +
    guides(color = guide_legend(override.aes = list(alpha = 3) ) )
ggsave('Figure 1.3.2b - family A - Adj_AF vs PHRED.png', plot = last_plot(), device = 'png', path = paste0(getwd(),plot.path), scale = 1, height = plot_height, width = plot_width, dpi = 'screen', units = 'in', limitsize = TRUE)    

#figure 1.3.3a
ggplot(only.family.B.dt[GNOMAD_AF!='.'],aes(x=LRFR, y=PHRED, color=Color)) +
    labs(x='Log Relative Frequency Ratio',y='PHRED') +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_point(alpha = 0.2,size = geom_size) + 
    scale_color_manual(name = "", breaks = c('Rare', 'Important','Impactful'), 
                       values=c('black','orange','red','gold'), guide=guide_legend(reverse=TRUE)) +
    guides(color = guide_legend(override.aes = list(alpha = 3) ) )
ggsave('Figure 1.3.3a - family B - LFRF vs PHRED.png', plot = last_plot(), device = 'png', path = paste0(getwd(),plot.path), scale = 1, height = plot_height, width = plot_width, dpi = 'screen', units = 'in', limitsize = TRUE)    

#figure 1.3.3b
ggplot(only.family.B.dt[GNOMAD_AF=='.'],aes(x=adj_AF, y=PHRED, color=Color)) +
    labs(x='Adjusted Cohort Allele Frequency',y='PHRED') +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_point(alpha = 0.2,size = geom_size) + 
    scale_color_manual(name = "", breaks = c('Rare', 'Important','Impactful'), 
                       values=c('gold','red'), guide=guide_legend(reverse=TRUE)) +
    guides(color = guide_legend(override.aes = list(alpha = 3) ) )
ggsave('Figure 1.3.3b - family B - Adj_AF vs PHRED.png', plot = last_plot(), device = 'png', path = paste0(getwd(),plot.path), scale = 1, height = plot_height, width = plot_width, dpi = 'screen', units = 'in', limitsize = TRUE)    

#create the tables
dir.create(paste0(getwd(),'/Results/germline'), showWarnings = FALSE)
write.csv(only.family.A.dt, paste0(getwd(),'/Results/germline/only.familyA.csv'), row.names = FALSE)
write.csv(only.family.B.dt, paste0(getwd(),'/Results/germline/only.familyB.csv'), row.names = FALSE)
write.csv(all.five.dt, paste0(getwd(),'/Results/germline/all.five.csv'), row.names = FALSE)
write.csv(final.w.CADD.cosmic, paste0(getwd(),'/Results/germline/final.w.CADD.cosmic.csv'), row.names = FALSE)
write.csv(pathogenic.cosmic.dt, paste0(getwd(),'/Results/germline/pathogenic.cosmic.csv'), row.names = FALSE)