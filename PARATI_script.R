# install packages if not
if(!require(dplyr)){install.packages("dplyr")}
if(!require(data.table)){install.packages("data.table")}
if(!require(stringr)){install.packages("stringr")}
if(!require(optparse)){install.packages("optparse")}
if(!require(openxlsx)){install.packages("openxlsx")}
if(!require(vcfR)){install.packages("vcfR")}

# Load necessary libraries
library(dplyr)
library(data.table)
library(stringr)
library(optparse)
library(openxlsx)
library(vcfR)


# read in the arguments that user input
option_list <- list(
  make_option(c("-g", "--geno"), type="character", help="Genotype vcf file in .vcf.gz", metavar="FILE"), # required
  make_option(c("-f", "--family"), type="character", help="Family ID Matching file in .xlsx", metavar="FILE"), # required
  make_option(c("-o", "--out"), type="character", help="Directory of output", metavar="PREFIX"), # required
  make_option(c("-c", "--chr"), type="integer", help="Chromosome", metavar="N"), # required
  make_option(c("-l","--haplen", type="integer", default=500000, help="Haplotype length for heterozygotes inference",metavar="N")), # optional
  make_option(c("-s","--savetemp", type="logical", default = FALSE, help="T or F, save temp directory or not, default is F")), # optional
  make_option(c("-b","--makebed", type="logical", default = FALSE, help="T or F, save in plink format or not, default is F")) # optional
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


req <- c("geno", "family", "out", "chr")
miss <- req[!nzchar(ifelse(sapply(req, function(x) as.character(opt[[x]])), "", ""))]
if (any(vapply(req, function(x) is.null(opt[[x]]), logical(1)))) {
  print_help(OptionParser(option_list=option_list))
  stop("Missing required option(s): ", paste(req[vapply(req, function(x) is.null(opt[[x]]), logical(1))], collapse=", "))
}


# Use the arguments
geno_file <- opt$geno
fam_file <- opt$family
dir_out <- opt$out 
chr_temp <- as.integer(opt$chr)


if (is.null(opt$haplen)) {hap_length=500000}
if (!is.null(opt$haplen)) {hap_length <- as.integer(opt$haplen)}
if (hap_length==500000) {print("Default Haplotype Length 500kb is used")}
if (!hap_length==500000) {print(paste0("Haplotype Length ", hap_length, " is used"))}

if (is.null(opt$savetemp)) {temp=F}
if (!is.null(opt$savetemp)) {temp <- as.logical(opt$savetemp)}
if (temp) {print(paste0("Temp Directory ", dir_out, "temp/ will be saved"))}

if (is.null(opt$makebed)) {plink_format=F}
if (!is.null(opt$makebed)) {plink_format=as.logical(opt$makebed)}
if (plink_format==T) {print("The output file will be saved in plink bed format")}
if (plink_format==F) {print("The output file will be saved in vcf format")}

print(hap_length)
print(temp)
print(plink_format)

###############################################################################################
dir_temp=paste0(dir_out,"temp/")
dir.create(dir_temp)
dir_vcf_trans=paste0(dir_temp,"vcf_Trios_ParentalTrans_WithHomozygote/")
system(paste0("mkdir -p ",dir_vcf_trans))
dir_vcf_nontrans=paste0(dir_temp,"vcf_Trios_ParentalNonTrans_WithHomozygote/")
system(paste0("mkdir -p ",dir_vcf_nontrans))
dir_vcf_HaplotypeMatching=paste0(dir_temp,"vcf_Trios_HaplotypeMatchingForHeterozygote/")
system(paste0("mkdir -p ",dir_vcf_HaplotypeMatching))
###############################################################################################



###############################################################################################
# identify parental transmitted vs non-transmitted allele by ethnicity and by chr
print(paste0("read in trios genotype info ", geno_file))
vcf_all=fread(file=paste0(geno_file))
setnames(vcf_all,"#CHROM","CHROM")
vcf_ByChr=vcf_all[CHROM==chr_temp]
setnames(vcf_ByChr,"CHROM","#CHROM")


info_column=names(vcf_ByChr)[1:9]
vcf_working=copy(vcf_ByChr)
setnames(vcf_working,old=names(vcf_working),new=str_split_fixed(names(vcf_working),"_",2)[,1])

print(paste0("read in family ID and role ", fam_file))
FID_linkerID=as.data.table(read.xlsx(paste0(fam_file)))
list_FID=unique(FID_linkerID$FamilyIndex)
print(paste0("There are ", length(list_FID), " families for inference"))


# extract vcf_info columns
vcf_info=vcf_working[,colnames(vcf_working)%in%c(info_column),with=F]
vcf_trans=copy(vcf_info)
vcf_nontrans=copy(vcf_info)
df_sim_perc_summary_AllFID=NA


for (FID_temp in list_FID) {
  print(FID_temp)
  fam_temp = FID_linkerID[FamilyIndex==FID_temp] |>  #get the original ID for a family of 3
    mutate(Role_BMP = case_when(Role=="C"~"B",
                                Role=="M"~"M",
                                Role=="F"~"P")) 
  iid = fam_temp$IndividualID
  member_id = fam_temp$Role_BMP   
  # examine haplotype
  vcf_temp=vcf_working[,colnames(vcf_working)%in%c(info_column)|colnames(vcf_working)%in%iid,with=F]
  setnames(vcf_temp,old=iid,new=member_id)  # change colnames of the trio to B M P
  
  vcf_temp[,B:=gsub("/","|",B)]
  vcf_temp[,M:=gsub("/","|",M)]
  vcf_temp[,P:=gsub("/","|",P)]
  
  vcf_temp[,B_hap1:=str_split_fixed(B, "[|]", 2)[,1]]
  vcf_temp[,B_hap2:=str_split_fixed(B, "[|]", 2)[,2]]
  vcf_temp[,M_hap1:=str_split_fixed(M, "[|]", 2)[,1]]
  vcf_temp[,M_hap2:=str_split_fixed(M, "[|]", 2)[,2]]
  vcf_temp[,P_hap1:=str_split_fixed(P, "[|]", 2)[,1]]
  vcf_temp[,P_hap2:=str_split_fixed(P, "[|]", 2)[,2]]
  
  unique(vcf_temp[,.(B_hap1,B_hap2,M_hap1,M_hap2,P_hap1,P_hap2)])
  unique(vcf_temp[,.(B_hap1,B_hap2)])
  unique(vcf_temp[,.(M_hap1,M_hap2)])
  unique(vcf_temp[,.(P_hap1,P_hap2)])
  unique(vcf_temp[,.(B)])
  unique(vcf_temp[,.(M)])
  unique(vcf_temp[,.(P)])
  
  # add columns for maternal and paternal transmitted vs non-transmitted alleles
  # update the values in these columns based on the combination of maternal, paternal, and child genotype/haplotype
  # based on the method described in PMID26284790
  vcf_temp[,M_transmitted:=NA]
  vcf_temp[,M_nontransmitted:=NA]
  vcf_temp[,P_transmitted:=NA]
  vcf_temp[,P_nontransmitted:=NA]
  
  ##################################################
  # Important assumption for homozygotes
  # in this script, we consider one of the alleles being transmitted and one of the alleles being non-transmitted for homozygotes
  # e.g., in case of "1|1", 1 is both transmitted and non-transmitted
  ##################################################
  
  # If mother, father, and child are homozygotes, the allele transmission can be unambiguously determined from their genotypes
  # check if there is any phasing error 
  vcf_temp[M=="1|1"&P=="1|1",]
  vcf_temp[M=="0|0"&P=="0|0",]
  unique(vcf_temp[M=="1|1"&P=="1|1",]$B)
  unique(vcf_temp[M=="0|0"&P=="0|0",]$B)
  
  # consider one of the alleles being transmitted for homozygotes
  vcf_temp[,M_transmitted:=ifelse(M=="1|1",1,M_transmitted)]
  vcf_temp[,M_transmitted:=ifelse(M=="0|0",0,M_transmitted)]
  vcf_temp[,P_transmitted:=ifelse(P=="1|1",1,P_transmitted)]
  vcf_temp[,P_transmitted:=ifelse(P=="0|0",0,P_transmitted)]
  
  # consider one of the alleles being non-transmitted for homozygotes (WithHomozygote)
  vcf_temp[,M_nontransmitted:=ifelse(M=="1|1","1",M_nontransmitted)]
  vcf_temp[,M_nontransmitted:=ifelse(M=="0|0","0",M_nontransmitted)]
  vcf_temp[,P_nontransmitted:=ifelse(P=="1|1","1",P_nontransmitted)]
  vcf_temp[,P_nontransmitted:=ifelse(P=="0|0","0",P_nontransmitted)]
  
  # If mother is homozygote but father is heterozygote, check the child's genotype
  # M_nontransmitted and M_transmitted have been updated in the previous step (maternal homozygote)
  # here, we determine paternal transmitted and non-transmitted allele 
  vcf_temp[(M=="1|1")&(P=="0|1"|P=="1|0"),] 
  unique(vcf_temp[(M=="1|1")&(P=="0|1"|P=="1|0"),][,.(B,M,P)])
  unique(vcf_temp[(M=="1|1")&(P=="0|1"|P=="1|0"),][,.(B_hap1,B_hap2,M_hap1,M_hap2,P_hap1,P_hap2)])
  
  vcf_temp[(M=="0|0")&(P=="0|1"|P=="1|0"),] 
  unique(vcf_temp[(M=="0|0")&(P=="0|1"|P=="1|0"),][,.(B,M,P)])
  unique(vcf_temp[(M=="0|0")&(P=="0|1"|P=="1|0"),][,.(B_hap1,B_hap2,M_hap1,M_hap2,P_hap1,P_hap2)])
  
  vcf_temp[,P_transmitted:=ifelse((M=="1|1")&(P=="0|1"|P=="1|0")&(B=="0|1"|B=="1|0"),0,P_transmitted)]
  vcf_temp[,P_transmitted:=ifelse((M=="1|1")&(P=="0|1"|P=="1|0")&(B=="1|1"),1,P_transmitted)]
  vcf_temp[,P_nontransmitted:=ifelse((M=="1|1")&(P=="0|1"|P=="1|0")&(B=="0|1"|B=="1|0"),1,P_nontransmitted)]
  vcf_temp[,P_nontransmitted:=ifelse((M=="1|1")&(P=="0|1"|P=="1|0")&(B=="1|1"),0,P_nontransmitted)]
  
  vcf_temp[,P_transmitted:=ifelse((M=="0|0")&(P=="0|1"|P=="1|0")&(B=="0|1"|B=="1|0"),1,P_transmitted)]
  vcf_temp[,P_transmitted:=ifelse((M=="0|0")&(P=="0|1"|P=="1|0")&(B=="0|0"),0,P_transmitted)]
  vcf_temp[,P_nontransmitted:=ifelse((M=="0|0")&(P=="0|1"|P=="1|0")&(B=="0|1"|B=="1|0"),0,P_nontransmitted)]
  vcf_temp[,P_nontransmitted:=ifelse((M=="0|0")&(P=="0|1"|P=="1|0")&(B=="0|0"),1,P_nontransmitted)]
  
  # If father is homozygote but mother is heterozygote, check the child's genotype
  # P_nontransmitted and P_transmitted have been updated in the previous step (paternal homozygote)
  #  here, we determine maternal transmitted and non-transmitted allele 
  vcf_temp[(M=="0|1"|M=="1|0")&(P=="1|1"),] 
  unique(vcf_temp[(M=="0|1"|M=="1|0")&(P=="1|1"),][,.(B,M,P)])
  unique(vcf_temp[(M=="0|1"|M=="1|0")&(P=="1|1"),][,.(B_hap1,B_hap2,M_hap1,M_hap2,P_hap1,P_hap2)])
  
  vcf_temp[(M=="0|1"|M=="1|0")&(P=="0|0"),] 
  unique(vcf_temp[(M=="0|1"|M=="1|0")&(P=="0|0"),][,.(B,M,P)])
  unique(vcf_temp[(M=="0|1"|M=="1|0")&(P=="0|0"),][,.(B_hap1,B_hap2,M_hap1,M_hap2,P_hap1,P_hap2)])
  
  vcf_temp[,M_transmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="1|1")&(B=="0|1"|B=="1|0"),0,M_transmitted)]
  vcf_temp[,M_transmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="1|1")&(B=="1|1"),1,M_transmitted)]
  vcf_temp[,M_nontransmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="1|1")&(B=="0|1"|B=="1|0"),1,M_nontransmitted)]
  vcf_temp[,M_nontransmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="1|1")&(B=="1|1"),0,M_nontransmitted)]
  
  vcf_temp[,M_transmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="0|0")&(B=="0|1"|B=="1|0"),1,M_transmitted)]
  vcf_temp[,M_transmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="0|0")&(B=="0|0"),0,M_transmitted)]
  vcf_temp[,M_nontransmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="0|0")&(B=="0|1"|B=="1|0"),0,M_nontransmitted)]
  vcf_temp[,M_nontransmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="0|0")&(B=="0|0"),1,M_nontransmitted)]
  
  # If mother and father are heterozygotes, and child is homozygote, the allele transmission can be unambiguously determined from their genotypes
  vcf_temp[(M=="0|1"|M=="1|0")&(P=="0|1"|P=="1|0")&(B=="1|1"|B=="0|0"),] 
  unique(vcf_temp[(M=="0|1"|M=="1|0")&(P=="0|1"|P=="1|0")&(B=="1|1"|B=="0|0"),][,.(B,M,P)])
  
  vcf_temp[,M_transmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="0|1"|P=="1|0")&(B=="1|1"),1,M_transmitted)]
  vcf_temp[,M_nontransmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="0|1"|P=="1|0")&(B=="1|1"),0,M_nontransmitted)]
  vcf_temp[,M_transmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="0|1"|P=="1|0")&(B=="0|0"),0,M_transmitted)]
  vcf_temp[,M_nontransmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="0|1"|P=="1|0")&(B=="0|0"),1,M_nontransmitted)]
  
  vcf_temp[,P_transmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="0|1"|P=="1|0")&(B=="1|1"),1,P_transmitted)]
  vcf_temp[,P_nontransmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="0|1"|P=="1|0")&(B=="1|1"),0,P_nontransmitted)]
  vcf_temp[,P_transmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="0|1"|P=="1|0")&(B=="0|0"),0,P_transmitted)]
  vcf_temp[,P_nontransmitted:=ifelse((M=="0|1"|M=="1|0")&(P=="0|1"|P=="1|0")&(B=="0|0"),1,P_nontransmitted)]
  
  
  
  ##################################################
  # If mother, father, and child are heterozygotes, we constructed a long-range (up to 1 Mb, +/-500kb) local haplotype (per Chr) around the SNP under consideration 
  # and compared haplotype sharing to determine allelic transmission (as described in PMID26284790)
  vcf_temp_hetero=vcf_temp[(M=="0|1"|M=="1|0")&(P=="0|1"|P=="1|0")&(B=="0|1"|B=="1|0"),] 
  vcf_temp_nonhetero=vcf_temp[!POS%in%vcf_temp_hetero$POS,]
  unique(vcf_temp_hetero[,.(B,M,P)])
  unique(vcf_temp_nonhetero[,.(B,M,P)])
  
  if (nrow(vcf_temp_hetero)>0){
    # Set up functions for analyses
    func.Parallel = function(i) {  
      # target variant
      vcf_temp_hetero_i=vcf_temp_hetero[i,]
      POS_temp=as.integer(unlist(vcf_temp_hetero_i[,.(POS)]))
      # extract haplotype within the +/- bpwindow from the target variant
      bpwindow=hap_length
      vcf_temp_POS1MB=vcf_temp[POS>POS_temp-bpwindow&POS<POS_temp+bpwindow,]
      # keep only the position with complete genotype data for the trios
      vcf_temp_POS1MB=vcf_temp_POS1MB[B_hap1%in%c(1,0)&B_hap2%in%c(1,0)&
                                        M_hap1%in%c(1,0)&M_hap2%in%c(1,0)&
                                        P_hap1%in%c(1,0)&P_hap2%in%c(1,0),]
      # identify the haplotype for each strand from the trios
      B_hap1_collapse <- paste(vcf_temp_POS1MB$B_hap1, collapse = "")
      B_hap2_collapse <- paste(vcf_temp_POS1MB$B_hap2, collapse = "")
      M_hap1_collapse <- paste(vcf_temp_POS1MB$M_hap1, collapse = "")
      M_hap2_collapse <- paste(vcf_temp_POS1MB$M_hap2, collapse = "")
      P_hap1_collapse <- paste(vcf_temp_POS1MB$P_hap1, collapse = "")
      P_hap2_collapse <- paste(vcf_temp_POS1MB$P_hap2, collapse = "")
      
      # check if the haplotype match between child and parents
      B_hap1_collapse %in% c(M_hap1_collapse, M_hap2_collapse)
      B_hap2_collapse %in% c(M_hap1_collapse, M_hap2_collapse)
      B_hap1_collapse %in% c(P_hap1_collapse, P_hap2_collapse)
      B_hap2_collapse %in% c(P_hap1_collapse, P_hap2_collapse)
      
      # usually, identical haplotypes between child and parents are rarely found
      # so we use a similarity percentage to compare each pair of haplotypes (child vs parents)
      # e.g., a high similarity percentage between B_hap1 and M_hap1 indicates that B_hap1 is more likely to be inherited from the mother (M_hap1)
      # thus for the target variant (POS_temp), we consider M_hap1 the transmitted allele and M_hap2 the non-transmitted allele 
      B_hap1_vs_M_hap1 <- mean(strsplit(B_hap1_collapse, "")[[1]] == strsplit(M_hap1_collapse, "")[[1]])
      B_hap1_vs_M_hap2 <- mean(strsplit(B_hap1_collapse, "")[[1]] == strsplit(M_hap2_collapse, "")[[1]])
      B_hap1_vs_P_hap1 <- mean(strsplit(B_hap1_collapse, "")[[1]] == strsplit(P_hap1_collapse, "")[[1]])
      B_hap1_vs_P_hap2 <- mean(strsplit(B_hap1_collapse, "")[[1]] == strsplit(P_hap2_collapse, "")[[1]])
      
      B_hap2_vs_M_hap1 <- mean(strsplit(B_hap2_collapse, "")[[1]] == strsplit(M_hap1_collapse, "")[[1]])
      B_hap2_vs_M_hap2 <- mean(strsplit(B_hap2_collapse, "")[[1]] == strsplit(M_hap2_collapse, "")[[1]])
      B_hap2_vs_P_hap1 <- mean(strsplit(B_hap2_collapse, "")[[1]] == strsplit(P_hap1_collapse, "")[[1]])
      B_hap2_vs_P_hap2 <- mean(strsplit(B_hap2_collapse, "")[[1]] == strsplit(P_hap2_collapse, "")[[1]])
      
      df_sim_perc=as.data.frame(t(data.frame(B_hap1_vs_M_hap1,B_hap1_vs_M_hap2,B_hap1_vs_P_hap1,B_hap1_vs_P_hap2,
                                             B_hap2_vs_M_hap1,B_hap2_vs_M_hap2,B_hap2_vs_P_hap1,B_hap2_vs_P_hap2)))
      df_sim_perc$pair=rownames(df_sim_perc)
      df_sim_perc=data.table(df_sim_perc)
      df_sim_perc[,B_hap:=str_split_fixed(pair,"_vs_",2)[,1]]
      df_sim_perc[,PM_hap:=str_split_fixed(pair,"_vs_",2)[,2]]
      setnames(df_sim_perc,old="V1",new="sim_perc")
      
      df_sim_perc[,nSNP_haplotype:=nrow(vcf_temp_POS1MB)]
      
      # if the highest sim_perc is smaller than 0.7, we consider low confidence for the allele transmission inference based on haplotype
      # and we assign NA for P_transmitted, P_nontransmitted, M_transmitted, M_nontransmitted for this target variant
      df_sim_perc=df_sim_perc[order(sim_perc,decreasing = TRUE),]
      if (max(df_sim_perc$sim_perc,na.rm = T)<0.7) {
        vcf_temp_hetero_i[,P_transmitted:=NA]
        vcf_temp_hetero_i[,P_nontransmitted:=NA]
        vcf_temp_hetero_i[,M_transmitted:=NA]
        vcf_temp_hetero_i[,M_nontransmitted:=NA]
        
        df_sim_perc[,status:="Low similarity"]
        
      } else {
        # determine the transimitted and non-transmitted alleles for the parent in the pair of the highest sim_perc
        df_sim_perc_max=df_sim_perc[sim_perc==max(sim_perc,na.rm = T),]
        PM_hap_temp=df_sim_perc_max$PM_hap
        
        # inference is only performed when there is a definite highest sim_perc 
        # i.e., when there are more than a single highest sim_perc, we consider this an ambiguous scenario
        if (length(PM_hap_temp)==1) {
          
          if (PM_hap_temp%like%"P_hap") {
            vcf_temp_hetero_i[,P_transmitted:=ifelse(PM_hap_temp=="P_hap1",P_hap1,
                                                     ifelse(PM_hap_temp=="P_hap2",P_hap2,NA))]
            vcf_temp_hetero_i[,P_nontransmitted:=ifelse(PM_hap_temp=="P_hap1",P_hap2,
                                                        ifelse(PM_hap_temp=="P_hap2",P_hap1,NA))]
          } else if (PM_hap_temp%like%"M_hap") {
            vcf_temp_hetero_i[,M_transmitted:=ifelse(PM_hap_temp=="M_hap1",M_hap1,
                                                     ifelse(PM_hap_temp=="M_hap2",M_hap2,NA))]
            vcf_temp_hetero_i[,M_nontransmitted:=ifelse(PM_hap_temp=="M_hap1",M_hap2,
                                                        ifelse(PM_hap_temp=="M_hap2",M_hap1,NA))]
          }
          
          # given that we determined the transimitted and non-transmitted alleles for one of the parents in the above step (based on the highest sim_perc)
          # the transimitted and non-transmitted alleles for the other parent are fixed (given that the trios are all heterozygote at this target variant)
          if (PM_hap_temp%like%"P_hap") {
            vcf_temp_hetero_i[,M_transmitted:=P_nontransmitted]
            vcf_temp_hetero_i[,M_nontransmitted:=P_transmitted]
          } else if (PM_hap_temp%like%"M_hap") {
            vcf_temp_hetero_i[,P_transmitted:=M_nontransmitted]
            vcf_temp_hetero_i[,P_nontransmitted:=M_transmitted]
          }
          
          df_sim_perc[,status:="Inferred based on haplotype"]
          
        } else if (length(PM_hap_temp)>1) {
          # we consider this an ambiguous scenario and no inference is performed
          vcf_temp_hetero_i[,P_transmitted:=NA]
          vcf_temp_hetero_i[,P_nontransmitted:=NA]
          vcf_temp_hetero_i[,M_transmitted:=NA]
          vcf_temp_hetero_i[,M_nontransmitted:=NA]
          
          df_sim_perc[,status:="Ambiguous"]
        }
        
      }
      
      df_sim_perc[,bpwindow:=bpwindow]
      df_sim_perc[,order:=c(1:nrow(df_sim_perc))]
      df_sim_perc=cbind(vcf_temp_hetero_i[,c("#CHROM","POS","ID"),with=F],df_sim_perc)
      
      vcf.summary <- data.frame(i, 
                                vcf_temp_hetero_i)
      invisible(vcf.summary)
      
      output <- list(df_sim_perc = df_sim_perc, vcf.summary = vcf.summary)
      return(output)
    }
    
    # run
    # vcf_temp_hetero_test=vcf_temp_hetero[1:100,]
    # List.result<-lapply(setNames(seq_len(nrow(vcf_temp_hetero_test)), c(1:nrow(vcf_temp_hetero_test))), func.Parallel)
    # system.time(List.result<-lapply(setNames(seq_len(nrow(vcf_temp_hetero_test)), c(1:nrow(vcf_temp_hetero_test))), func.Parallel))
    List.result<-lapply(setNames(seq_len(nrow(vcf_temp_hetero)), c(1:nrow(vcf_temp_hetero))), 
                        func.Parallel)
    
    # export similarity percentage for all target haplotype
    df_sim_perc_summary <- rbindlist(lapply(List.result, function(x) x[[1]]))
    table(unique(df_sim_perc_summary[,.(ID,status)])$status)
    df_sim_perc_summary[order==1,][order(sim_perc),]
    df_sim_perc_summary[order==2,][order(sim_perc),]
    df_sim_perc_summary[,FID:=gsub("_","-",FID_temp)]
    
    # summarise transmitted and non-transmitted allele for the heterozygote
    vcf_temp_hetero_summary <- rbindlist(lapply(List.result, function(x) x[[2]]))
    setnames(vcf_temp_hetero_summary,old="X.CHROM",new="#CHROM")
    vcf_temp_hetero_summary[,i:=NULL]
    
  }
  ##################################################
  
  
  
  ##################################################
  # combine allele transmission from heterozygote inference and from non-heterozygote
  vcf_temp=rbind(vcf_temp_hetero_summary,vcf_temp_nonhetero)
  
  # examine maternal and paternal transmitted vs non-transmitted alleles not yet updated
  # the remaining NAs are likely due to unmatched hyplotype in the +/- bpwindow
  # matching is possible with a narrower bpwindow 
  vcf_temp[is.na(M_transmitted)==T|is.na(M_nontransmitted)==T|
             is.na(P_transmitted)==T|is.na(P_nontransmitted)==T,]
  
  unique(vcf_temp[is.na(M_transmitted)==T|is.na(M_nontransmitted)==T|
                    is.na(P_transmitted)==T|is.na(P_nontransmitted)==T,][,.(B,M,P)])
  
  vcf_temp[,M_transmitted:=ifelse(is.na(M_transmitted)==T,".",M_transmitted)]
  vcf_temp[,M_nontransmitted:=ifelse(is.na(M_nontransmitted)==T,".",M_nontransmitted)]
  vcf_temp[,P_transmitted:=ifelse(is.na(P_transmitted)==T,".",P_transmitted)]
  vcf_temp[,P_nontransmitted:=ifelse(is.na(P_nontransmitted)==T,".",P_nontransmitted)]
  
  # re-generate format for vcf.gz, i.e., "x|x" or "x|." format
  vcf_temp[,M_transmitted_vcf:=paste0(M_transmitted,"|",".")]
  vcf_temp[,M_nontransmitted_vcf:=paste0(M_nontransmitted,"|",".")]
  vcf_temp[,P_transmitted_vcf:=paste0(P_transmitted,"|",".")]
  vcf_temp[,P_nontransmitted_vcf:=paste0(P_nontransmitted,"|",".")]
  
  # extract columns for corresponding outputs
  vcf_temp_trans=vcf_temp[,colnames(vcf_temp)%in%c(info_column)|colnames(vcf_temp)%like%"_transmitted_vcf",with=F]
  setnames(vcf_temp_trans,
           old=c("M_transmitted_vcf","P_transmitted_vcf"),
           #new=c(paste0(c("M","P"),gsub("_","-",FID_temp)))
           new=c(fam_temp[Role=="M"]$IndividualID, fam_temp[Role=="F"]$IndividualID))
  vcf_trans=merge(vcf_trans,vcf_temp_trans,by=info_column,all.x=T)
  
  vcf_temp_nontrans=vcf_temp[,colnames(vcf_temp)%in%c(info_column)|colnames(vcf_temp)%like%"_nontransmitted_vcf",with=F]
  setnames(vcf_temp_nontrans,
           old=c("M_nontransmitted_vcf","P_nontransmitted_vcf"),
           #new=c(paste0(c("M","P"),gsub("_","-",FID_temp)))
           new=c(fam_temp[Role=="M"]$IndividualID, fam_temp[Role=="F"]$IndividualID))
  vcf_nontrans=merge(vcf_nontrans,vcf_temp_nontrans,by=info_column,all.x=T)
  
  # combine df_sim_perc_summary for all FID
  df_sim_perc_summary_AllFID=rbind(df_sim_perc_summary_AllFID,df_sim_perc_summary,fill=T)
}


# save the inferred genotype
write.table(vcf_trans,file=gzfile(paste0(dir_vcf_trans,"trans_genotype_chr_",chr_temp,".vcf.gz")),row.names = F,quote=F,sep="\t")
write.table(vcf_nontrans,file=gzfile(paste0(dir_vcf_nontrans,"nontrans_genotype_chr_",chr_temp,".vcf.gz")),row.names = F,quote=F,sep="\t")

# save the heterozygote matching reference
df_sim_perc_summary_AllFID[,x:=NULL]
df_sim_perc_summary_AllFID=df_sim_perc_summary_AllFID[is.na(FID)==F,]
df_sim_perc_summary_AllFID=df_sim_perc_summary_AllFID[,c("FID","#CHROM","POS","ID","pair","B_hap","PM_hap","bpwindow","nSNP_haplotype","sim_perc","order","status")]
write.csv(df_sim_perc_summary_AllFID,file=gzfile(paste0(dir_vcf_HaplotypeMatching,"/Trios_chr",chr_temp,"_HaplotypeMatchingForHeterozygote.csv.gz")),row.names = F)






########################################################################################
# change the format to real vcf format
df_nontrans <- fread(paste0(dir_vcf_nontrans,"nontrans_genotype_chr_",chr_temp,".vcf.gz"))
df_trans <- fread(paste0(dir_vcf_trans,"trans_genotype_chr_",chr_temp,".vcf.gz"))

vcf_full <- read.vcfR(paste0(geno_file))
vcf_full_chr <- vcf_full[ getCHROM(vcf_full) %in% c(as.character(chr_temp),paste0("chr",chr_temp)), ]
meta_section <- vcf_full_chr@meta


vcf_nontrans_eth <- new("vcfR",
                        meta = meta_section,
                        fix = as.matrix(df_nontrans[, 1:8]),
                        gt = as.matrix(df_nontrans[,9:ncol(df_nontrans)]))
write.vcf(vcf_nontrans_eth, file=paste0(dir_vcf_nontrans,"chr",chr_temp,".vcf"))

vcf_trans_eth <- new("vcfR",
                      meta = meta_section,
                      fix = as.matrix(df_trans[, 1:8]),
                      gt = as.matrix(df_trans[,9:ncol(df_trans)]))
write.vcf(vcf_trans_eth, file=paste0(dir_vcf_trans,"chr",chr_temp,".vcf"))




############################################################################
# change vcf to plink format
if (plink_format==T) {
  #system(paste0("module load plink/1.9"))
  vcf_nta=paste0(dir_vcf_nontrans,"chr",chr_temp,".vcf")
  out_nta=paste0(dir_out,"nontrans_chr_",chr_temp)
  system(paste0("plink --vcf ", vcf_nta,
                     " --make-bed ",
                     " --out ", out_nta,
                     " --vcf-half-call h"))
  
  vcf_ta=paste0(dir_vcf_trans,"chr",chr_temp,".vcf") 
  out_ta=paste0(dir_out,"trans_chr_",chr_temp)
  system(paste0("plink --vcf ", vcf_ta,
                " --make-bed ",
                " --out ", out_ta,
                " --vcf-half-call h"))
}

if (plink_format==F) {
  write.vcf(vcf_nontrans_eth, file=paste0(dir_out,"chr",chr_temp,".vcf")) # save in vcf format
  write.vcf(vcf_trans_eth, file=paste0(dir_out,"chr",chr_temp,".vcf"))
}


if (temp==F) {unlink(paste0(dir_temp), recursive=TRUE)}






