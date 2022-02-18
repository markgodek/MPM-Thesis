library(data.table)
library(stringr)

setwd('C:/Users/Mark/Dropbox (Partners HealthCare)/Mark Shared MPM Folder') #home 
#setwd('C:/Users/mi313/Dropbox (Partners HealthCare)/Mark Shared MPM Folder') #work

somatic.analysis <- function(germline.table, somatic.file){
  
  somatic.table<-process.somatic(somatic.file)
  
  keycols <- c('GeneID', 'patient')
  setkeyv(germline.table,keycols)
  setkeyv(somatic.table,keycols)
  return_table<-merge(x = somatic.table, y = germline.table, by = keycols, all.x = FALSE)
  
  colnames(return_table)[3]<-"Uploaded_variation"
  
  #create return table with germline columns on the lefthand side, sorted by CHROM and POS
  return_table<-return_table[
      with(return_table,order(as.numeric(CHROM),as.numeric(POS))),c(13:30,1:12)]
  
  #collapsing on GeneIDs is fine AFTER we intersect that data tables
  return_table[,GeneID:=paste(unique(GeneID), sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT", "Uploaded_variation", "Allele")]
  return_table[,Feature:=paste(unique(Feature), sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT", "Uploaded_variation", "Allele")]
  return_table[,cDNA_position:=paste(unique(cDNA_position), sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT", "Uploaded_variation", "Allele")]
  return_table[,Consequence :=paste(unique(Consequence) , sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT", "Uploaded_variation", "Allele")]
  return_table[,CDS_position:=paste(unique(CDS_position), sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT", "Uploaded_variation", "Allele")]
  return_table[,Protein_position:=paste(unique(Protein_position), sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT", "Uploaded_variation", "Allele")]
  return_table[,Amino_acids:=paste(unique(Amino_acids), sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT", "Uploaded_variation", "Allele")]
  return_table[,Codons:=paste(unique(Codons), sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT", "Uploaded_variation", "Allele")]
  return_table[,Location:=paste(unique(Location), sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT", "Uploaded_variation", "Allele")]

  return(unique(return_table))
}

process.somatic<-function(somatic.file){
  somatic.table<-as.data.table(fread(somatic.file, skip='#Uploaded_variation'))
  patient<-str_extract(somatic.file,'GNT[0-9][0-9][0-9]')
  colnames(somatic.table)[4]<-'GeneID'
  somatic.table$patient<-patient

  #somatic.table<-as.data.table(read.table(readLines(somatic.file), stringsAsFactors = FALSE))
  #print(head(somatic.table))
  return(unique(somatic.table[,!c('Feature_type','Extra','Existing_variation')]))

}

process.germline<-function(germline.file){
  
  AD_cutoff <- 15 #sets cutoff for insufficient allele depth
  DP_Cutoff <- 15 #sets cutoff for insufficient read depth 
  
  #read the germline alleles and melt the wide format into long format based on the columns listed
  af.dt<-fread(file=germline.file) 
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
  final.dt<-final.frame[,c('CHROM', 'POS', 'REF', 'ALT', 'AF', 'GNOMAD_AF', 'GENE', 'FEATURE', 'DP', 'AC', 'AN', 'sample','family.A.GT',
                           'family.B.GT','genotype_final')]
  final.dt<-final.dt[genotype_final!='./.'&genotype_final!='0/0'&genotype_final!='i/i'&genotype_final!='0/i']
  
  #store old allele frequency and calculate a new one based on the updated AC and AN
  final.dt[,AF_old:=AF]
  final.dt[,AF:=AC/AN]
  final.dt<-final.dt[AC>0&AN>0,]
  
  #add the family genotype data to the final data table
  setkeyv(family.frame,keycols)
  setkeyv(final.dt,keycols)
  final.dt<-merge(x=final.dt,y=family.frame,by=keycols,all.x = TRUE)
  
  
  #annotate samples with GeneIds
  CADD.dt<-fread(file = 'CADD_annotations_for_germline.PASS.tsv')
  colnames(CADD.dt)[1]<-'CHROM'
  colnames(CADD.dt)[2]<-'POS'
  colnames(CADD.dt)[3]<-'REF'
  colnames(CADD.dt)[4]<-'ALT'
  #CADD.dt[,GeneID:=paste(GeneID, sep = "", collapse = ", "),by=c("CHROM","POS","REF","ALT")]
  CADD.dt<-unique(CADD.dt[!is.na(GeneID)][,c('CHROM','POS','REF','ALT','GeneID','PHRED')])
  setkeyv(final.dt,keycols)
  setkeyv(CADD.dt,keycols)
  
  samples.w.geneIDs.germline<-merge(x = final.dt, y = CADD.dt, by = keycols, all.x = TRUE)
  samples.w.geneIDs.germline<-samples.w.geneIDs.germline[grep("ENSG[:digit:]*",GeneID)]  #remove rows without GeneIDs
  
  remove(family.frame,final.dt,CADD.dt,final.frame)
  
  #remove per-sample family genotype and rename per-variant family genotype column
  samples.w.geneIDs.germline<-samples.w.geneIDs.germline[,!c('family.A.GT.x','family.B.GT.x')]
  names(samples.w.geneIDs.germline)[names(samples.w.geneIDs.germline) == "family.A.GT.y"] <- "family.A.GT"
  names(samples.w.geneIDs.germline)[names(samples.w.geneIDs.germline) == "family.B.GT.y"] <- "family.B.GT"
  
  #add columns for genotypes of families
  samples.w.geneIDs.germline[,"A.0/0":=str_count(family.A.GT,'0/0')]
  samples.w.geneIDs.germline[,"A.0/i":=str_count(family.A.GT,'0/i')]
  samples.w.geneIDs.germline[,"A.0/1":=str_count(family.A.GT,'0/1')]
  samples.w.geneIDs.germline[,"A.1/1":=str_count(family.A.GT,'1/1')]
  samples.w.geneIDs.germline[,"A.i/1":=str_count(family.A.GT,'i/1')]
  samples.w.geneIDs.germline[,"A.i/i":=str_count(family.A.GT,'i/i')]
  samples.w.geneIDs.germline[,"A../.":=str_count(family.A.GT,'\\./\\.')]
  samples.w.geneIDs.germline[,"B.0/0":=str_count(family.B.GT,'0/0')]
  samples.w.geneIDs.germline[,"B.0/1":=str_count(family.B.GT,'0/1')]
  samples.w.geneIDs.germline[,"B.0/i":=str_count(family.B.GT,'0/i')]
  samples.w.geneIDs.germline[,"B.1/1":=str_count(family.B.GT,'1/1')]
  samples.w.geneIDs.germline[,"B.i/1":=str_count(family.B.GT,'i/1')]
  samples.w.geneIDs.germline[,"B.i/i":=str_count(family.B.GT,'i/i')]
  samples.w.geneIDs.germline[,"B../.":=str_count(family.B.GT,'\\./\\.')]
  
  #count number of zero or 1. if there are 2 or more i's, reduce observations to 1, else either two 1's or 0/1 were seen, i.e. 2 observations
  #add these observations from AN
  samples.w.geneIDs.germline[,'AN+.A':=ifelse(`A.1/1`>=1|`A.0/1`>=1|`A.0/0`>=1,2,ifelse(`A.i/1`>=1|`A.0/1`>=1,1,0))]
  samples.w.geneIDs.germline[,'AN+.B':=ifelse(`B.1/1`>=1|`B.0/1`>=1|`B.0/0`>=1,2,ifelse(`B.i/1`>=1|`B.0/1`>=1,1,0))]
  
  #subtract zero and 1 observations from AN and 1's from AC
  #if there is a genotype, count zeroes and 1's, else set it to zero - code to deal with NA's in genotype column
  samples.w.geneIDs.germline[,'0.counts':=str_count(family.A.GT,'0')+str_count(family.B.GT,'0')]
  samples.w.geneIDs.germline[,'1.counts':=str_count(family.A.GT,'1')+str_count(family.B.GT,'1')]
  
  #count number of 1's in the genotypes. 1/1 = 2 and overides 0/1 and i/1, which = 1
  samples.w.geneIDs.germline[,'AC+.A':=ifelse(`A.1/1`>=1,2,ifelse(`A.0/1`>=1|`A.i/1`>=1,1,0))]
  samples.w.geneIDs.germline[,'AC+.B':=ifelse(`B.1/1`>=1,2,ifelse(`B.0/1`>=1|`B.i/1`>=1,1,0))]
  
  #adjusted AN - AN with observations in all families flattened
  samples.w.geneIDs.germline[,'adj_AN':=AN+`AN+.A`+`AN+.B`-`1.counts`-`0.counts`]
  
  #adjusted AC - AC with all observations in families flattened
  samples.w.geneIDs.germline[,'adj_AC':=AC-`1.counts`+`AC+.A`+`AC+.B`]
  
  #adjusted AF - AF of variant after reducing observations per family to one
  samples.w.geneIDs.germline[,adj_AF:= adj_AC/adj_AN]
  
  #add some columns useful for filtering - high AF_distance means 
  samples.w.geneIDs.germline[,'family':=`A.1/1`+`A.0/1`+`A.i/1`+`B.1/1`+`B.0/1`+`B.i/1`]
  samples.w.geneIDs.germline[,'AF_dist':=abs(adj_AF-as.numeric(GNOMAD_AF))]
  samples.w.geneIDs.germline[,'LRFR':=log10(adj_AF/as.numeric(GNOMAD_AF))]
  samples.w.geneIDs.germline[,'Rare':=ifelse(LRFR>0.75|GNOMAD_AF=='.','Yes','No')]       # a variant is "rare" if it has a high log ratio, or if it does not appear in the gnomAD dataset
  samples.w.geneIDs.germline[,'Impactful':=ifelse(is.na(PHRED),'No',ifelse(PHRED>=20,'Yes','No'))]
  samples.w.geneIDs.germline[,'Important':=ifelse((Rare=='Yes'&Impactful=='Yes'),'Yes','No')]
  samples.w.geneIDs.germline[(GNOMAD_AF=='.'&PHRED>=20&adj_AF>0.125), 'Important':='Yes']         # Two samples that are important do not have GNOMAD_AF values, so annotated them manually
  
  #create a column to color the points in the plots with the heirarch important > impactful > rare
  samples.w.geneIDs.germline[,Classification:=ifelse(Important=='Yes','Important',ifelse(Impactful=='Yes','Impactful',ifelse(Rare=='Yes','Rare','Common')))]
  
  patients.w.geneIDs.germline <- patient.assigner(samples.w.geneIDs.germline)

  remove(A.frame, B.frame, family.A, family.B, samples.w.geneIDs.germline)
  
  returnCols<-c('CHROM','POS','REF','ALT','AF','AF_old','GNOMAD_AF','LRFR','DP','AC','AN',
                'GENE','FEATURE','sample','patient','genotype_final','GeneID','adj_AF','PHRED','Classification')
  
  return(patients.w.geneIDs.germline[,..returnCols])
}

patient.assigner <- function(table){
  
  #connect sample ID with patient ID
  GNT901<-c('HITS622838')
  GNT902<-c('HITS622839')
  GNT903<-c('HITS622840')
  GNT904<-c('HITS622841')
  GNT905<-c('HITS622842')
  GNT906<-c('HITS622843')
  GNT907<-c('HITS622844')
  GNT908<-c('HITS622845')
  GNT909<-c('HITS622847', 'HITS622846')
  GNT910<-c('HITS622849', 'HITS622848')
  GNT911<-c('HITS622851', 'HITS622850')
  GNT912<-c('HITS622853', 'HITS622852')
  GNT913<-c('HITS622855', 'HITS622854')
  GNT914<-c('HITS622856')
  GNT915<-c('HITS622858', 'HITS622857')
  GNT916<-c('HITS622860', 'HITS622859')
  GNT917<-c('HITS622862', 'HITS622861')
  GNT918<-c('HITS622864', 'HITS622863')
  GNT919<-c('HITS622866', 'HITS622865')
  GNT920<-c('HITS622868', 'HITS622867')
  GNT921<-c('HITS622870', 'HITS622869')
  GNT922<-c('HITS622872', 'HITS622871')
  GNT923<-c('HITS622875', 'HITS622873')
  GNT924<-c('HITS622877', 'HITS622876')
  
  table$patient[table$sample %in% GNT901] <- 'GNT901'
  table$patient[table$sample %in% GNT902] <- 'GNT902'
  table$patient[table$sample %in% GNT903] <- 'GNT903'
  table$patient[table$sample %in% GNT904] <- 'GNT904'
  table$patient[table$sample %in% GNT905] <- 'GNT905'
  table$patient[table$sample %in% GNT906] <- 'GNT906'
  table$patient[table$sample %in% GNT907] <- 'GNT907'
  table$patient[table$sample %in% GNT908] <- 'GNT908'
  table$patient[table$sample %in% GNT909] <- 'GNT909'
  table$patient[table$sample %in% GNT910] <- 'GNT910'
  table$patient[table$sample %in% GNT911] <- 'GNT911'
  table$patient[table$sample %in% GNT912] <- 'GNT912'
  table$patient[table$sample %in% GNT913] <- 'GNT913'
  table$patient[table$sample %in% GNT914] <- 'GNT914'
  table$patient[table$sample %in% GNT915] <- 'GNT915'
  table$patient[table$sample %in% GNT916] <- 'GNT916'
  table$patient[table$sample %in% GNT917] <- 'GNT917'
  table$patient[table$sample %in% GNT918] <- 'GNT918'
  table$patient[table$sample %in% GNT919] <- 'GNT919'
  table$patient[table$sample %in% GNT920] <- 'GNT920'
  table$patient[table$sample %in% GNT921] <- 'GNT921'
  table$patient[table$sample %in% GNT922] <- 'GNT922'
  table$patient[table$sample %in% GNT923] <- 'GNT923'
  table$patient[table$sample %in% GNT924] <- 'GNT924'
  
  return(table)
}

#germline.file<-"AF/germline/germline.PASS.AFonly.frq"            # for debugging process.germline function
germline.table<-process.germline("AF/germline/germline.PASS.AFonly.frq")

GNT909.overlap<-somatic.analysis(germline.table,"Results/Somatic/VEP Annotated/GNT909.geneID.PASS.txt")
GNT910.overlap<-somatic.analysis(germline.table,"Results/Somatic/VEP Annotated/GNT910.geneID.PASS.txt")
GNT911.overlap<-somatic.analysis(germline.table,"Results/Somatic/VEP Annotated/GNT911.geneID.PASS.txt")
GNT912.overlap<-somatic.analysis(germline.table,"Results/Somatic/VEP Annotated/GNT912.geneID.PASS.txt")
GNT913.overlap<-somatic.analysis(germline.table,"Results/Somatic/VEP Annotated/GNT913.geneID.PASS.txt")
GNT915.overlap<-somatic.analysis(germline.table,"Results/Somatic/VEP Annotated/GNT915.geneID.PASS.txt")
GNT916.overlap<-somatic.analysis(germline.table,"Results/Somatic/VEP Annotated/GNT916.geneID.PASS.txt")
GNT917.overlap<-somatic.analysis(germline.table,"Results/Somatic/VEP Annotated/GNT917.geneID.PASS.txt")
GNT918.overlap<-somatic.analysis(germline.table,"Results/Somatic/VEP Annotated/GNT918.geneID.PASS.txt")
GNT919.overlap<-somatic.analysis(germline.table,"Results/Somatic/VEP Annotated/GNT919.geneID.PASS.txt")
GNT920.overlap<-somatic.analysis(germline.table,"Results/Somatic/VEP Annotated/GNT920.geneID.PASS.txt")
GNT921.overlap<-somatic.analysis(germline.table,"Results/Somatic/VEP Annotated/GNT921.geneID.PASS.txt")
GNT922.overlap<-somatic.analysis(germline.table,"Results/Somatic/VEP Annotated/GNT922.geneID.PASS.txt")
GNT923.overlap<-somatic.analysis(germline.table,"Results/Somatic/VEP Annotated/GNT923.geneID.PASS.txt")
GNT924.overlap<-somatic.analysis(germline.table,"Results/Somatic/VEP Annotated/GNT924.geneID.PASS.txt")

germline.variants.with.somatic.hits.table<-do.call("rbind", list(GNT909.overlap,GNT910.overlap,GNT911.overlap,GNT912.overlap,
                                                                 GNT913.overlap,GNT915.overlap,GNT916.overlap,GNT917.overlap,
                                                                 GNT918.overlap,GNT919.overlap,GNT920.overlap,GNT921.overlap,
                                                                 GNT922.overlap,GNT923.overlap,GNT924.overlap))

write.csv(germline.variants.with.somatic.hits.table, paste0(getwd(),'/Results/Somatic/Germ-Som Overlap/germline.variants.with.somatic.hits.csv'), row.names = FALSE)


