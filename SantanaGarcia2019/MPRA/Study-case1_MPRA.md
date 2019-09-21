---
title: "Study-case1_MPRA"
author: "Maria Jose Rocha, Walter Santana"
date: "21/09/2019"
output: 
  html_document:
    keep_md: yes
---

## Description

In order to assess the quality of results retrieves by *variation-scan* we compared the results of our tool to the variant effect meassured by 
- Ulirsch, Jacob C., Satish K. Nandakumar, Li Wang, Felix C. Giani, Xiaolan Zhang, Peter Rogov, Alexandre Melnikov, et al. 2016. “Systematic Functional Dissection of Common Genetic Variation Affecting Red Blood Cell Traits.” Cell 165 (6): 1530–45.

Furthermore, we compared our results to the one obteined by the authors when they meassured the putative effect of the variants with the tools Deepsea and DeltaSVM.

## Dependecies
Several dependencies were used for data processing, in addition to R libraries. This tools
should also be present in the environmental variables.

* R (3.3.2)
* clustalw2 (2.1)
* htslib (1.9)
* samtools (1.9)
* bcftools (1.9)

## R Libraries
Libraries used throughout the code sections.

```r
# Biostrings(2.42.1)
library("Biostrings")
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter,
##     Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unsplit, which,
##     which.max, which.min
```

```
## Loading required package: S4Vectors
```

```
## Warning: package 'S4Vectors' was built under R version 3.6.1
```

```
## Loading required package: stats4
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## Loading required package: IRanges
```

```
## Warning: package 'IRanges' was built under R version 3.6.1
```

```
## Loading required package: XVector
```

```
## 
## Attaching package: 'Biostrings'
```

```
## The following object is masked from 'package:base':
## 
##     strsplit
```

```r
# ggplot2(3.0.0)
library("ggplot2")
library("readr")
library("preprocessCore")
library("reshape")
```

```
## 
## Attaching package: 'reshape'
```

```
## The following objects are masked from 'package:S4Vectors':
## 
##     expand, rename
```

```r
library("plyr")
```

```
## 
## Attaching package: 'plyr'
```

```
## The following objects are masked from 'package:reshape':
## 
##     rename, round_any
```

```
## The following object is masked from 'package:XVector':
## 
##     compact
```

```
## The following object is masked from 'package:IRanges':
## 
##     desc
```

```
## The following object is masked from 'package:S4Vectors':
## 
##     rename
```

```r
library("betareg")
library("qvalue")
library("UpSetR")
library("stringr")
library("ggpubr")
```

```
## Loading required package: magrittr
```

```
## 
## Attaching package: 'ggpubr'
```

```
## The following object is masked from 'package:plyr':
## 
##     mutate
```

```r
library("gclus")
```

```
## Loading required package: cluster
```

```r
library("pROC")
```

```
## Type 'citation("pROC")' for a citation.
```

```
## 
## Attaching package: 'pROC'
```

```
## The following objects are masked from 'package:IRanges':
## 
##     cov, var
```

```
## The following objects are masked from 'package:S4Vectors':
## 
##     cov, var
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     var
```

```
## The following objects are masked from 'package:stats':
## 
##     cov, smooth, var
```

## Declare used functions
R functions definitions.

```r
#Create consensus sequences for consMatrix 
getConsensusSeq <- function(consMatrix){
  consensusSequence <- apply(consMatrix, 2, function(column){
    select_cols <- column[column>0]
    select_cols <- select_cols[grepl("[^-]",names(select_cols))]
    if(length(select_cols) != 1){
      print( paste("Unable to proceed: a mismatch was found!",aligned_seq, collapse = "") )
      return(0)
    }
    return(names(select_cols))
  } )
  consensusSequence <- paste0(consensusSequence, collapse = "")
  return(consensusSequence)
}

#Test if sequences are continuous
check_seqcontinuity <- function(sequences) {
  for(aligned_seq in sequences) {
  countBases = 0
  countGaps  = 0
  for(base in strsplit(aligned_seq,"")[[1]]) {
    if ( base == "A" ) {
      countBases = countBases + 1
    } else if( base == "T" ) {
      countBases = countBases + 1
    } else if( base == "C"  ) {
      countBases = countBases + 1
    } else if( base == "G"  ) {
      countBases = countBases + 1
    } else if( base == "-"  ) {
      #If a gap is found between the sequences,it breaks
      if(countBases != 0 && countBases != 145) {
        print( paste("Unable to proceed: a gap was found between the aligned sequence ",aligned_seq, collapse = "") )
        return(0)
      }
      countGaps = countGaps + 1
    } else {
      print( paste("Unable to proceed: found base is not a A,T,C,G or -  in ",aligned_seq, collapse = "") )
      return(0)
    }
    
  }

  }
  return(1)
}
```

## Raw data download  

### MPRA assessed variants (Table S1)
The downloaded Table S1 file (TableS1_Constructs.csv), was obtained from https://www.cell.com/fulltext/S0092-8674(16)30493-7
(Ulirsch et al., 2016) in 15-06-2018 within the Supplemental Information section and then exported
as a CSV file. 

NOTE: It is important to pinpoint that in line 978 there was an extra column, i.e. 45939722	#REF!,
so we erased it before further processing. 


```bash
#Convert table to .tab
echo -e "#Variant\tConstructPos\tConstruct" > data/TableS1_Constructs.tab
perl -pe 's/\,/\t/g' data/TableS1_Constructs.csv | tail -n 16554 >> data/TableS1_Constructs.tab
#Convert table to vcf
echo -e "##fileformat=VCFv4.1" > data/MPRA_Variants.vcf
echo -e "##fileDate=20180615" >> data/MPRA_Variants.vcf
echo -e "##source=TableS1;url=https://www.cell.com/cell/fulltext/S0092-8674(16)30493-7" >> data/MPRA_Variants.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> data/MPRA_Variants.vcf
grep -v "#" data/TableS1_Constructs.tab | cut -f 1 | sort | uniq |awk \
'{split($1,variant,":"); print variant[1]"\t"variant[2]"\t.\t"variant[3]"\t"variant[4]"\t.\t.\t."}'| \
sort -k1,1 -k2,2n >> data/MPRA_Variants.vcf
```


### Ensembl release 75 variants and fasta genome
**Time consumming step, final files are provided**
The downloaded files from Ensembl were retrieved in 15-06-2018.

```bash
#Download ensembl variations release 75
wget --output-document=data/Homo_sapiens.r75.vcf.gz ftp://ftp.ensembl.org/pub/release-75/variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz
#Download release-75 fastas and concatenate them
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT
do
  wget --output-document=data/Homo_sapiens.GRCh37.75.dna.chromosome.$i.fa.gz ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.$i.fa.gz
  zcat < data/Homo_sapiens.GRCh37.75.dna.chromosome.$i.fa.gz|perl -pe s"/^>(.+)/>$i/" > data/Homo_sapiens.GRCh37.r75.chr$i.fa
  rm data/Homo_sapiens.GRCh37.75.dna.chromosome.$i.fa.gz
done
```

## Data Pre-processing

### MPRA variant annotation with Ensembl release 75 variations
The assessed MPRA variants were annotated using Ensembl release 75 variants. This entailed the retrieval of
the chromsome sizes from the GRCh37 human genome build used to map these variants. 

```bash
#Append all fasta files in a single file
cat data/Homo_sapiens.GRCh37.r75.chr*.fa| bgzip -c > data/Homo_sapiens.GRCh37.r75.fa.gz
rm data/Homo_sapiens.GRCh37.r75.chr*.fa

#Decompress
gzip -d data/Homo_sapiens.r75.vcf.gz 

#Compress to bgzip
bgzip -c data/Homo_sapiens.r75.vcf  > data/Homo_sapiens.r75.vcf.gz 
bgzip -c data/MPRA_Variants.vcf     > data/MPRA_Variants.vcf.gz

#Get fasta index
samtools faidx data/Homo_sapiens.GRCh37.r75.fa.gz

#Parse chromosome information and create a file with it
awk '{print "##contig=<ID="$1",length="$2">"}' data/Homo_sapiens.GRCh37.r75.fa.gz.fai > data/Homo_sapien.GRCh37.r75.txt

#Add chromosome information to vcfs
bcftools annotate --header-lines data/Homo_sapien.GRCh37.r75.txt --output data/Homo_sapiens.info.r75.vcf.gz --output-type z data/Homo_sapiens.r75.vcf.gz
bcftools annotate --header-lines data/Homo_sapien.GRCh37.r75.txt --output data/MPRA.info.vcf.gz --output-type z data/MPRA_Variants.vcf.gz

#Sort files
bcftools sort --output-file data/Homo_sapiens.r75.sort.vcf.gz --output-type z data/Homo_sapiens.info.r75.vcf.gz
bcftools sort --output-file data/MPRA_Variants.sort.vcf.gz    --output-type z data/MPRA.info.vcf.gz

#Index files
tabix -p vcf data/Homo_sapiens.r75.sort.vcf.gz
tabix -p vcf data/MPRA_Variants.sort.vcf.gz

#Annotate variants
bcftools annotate --output data/MPRA_Variants.annot75.vcf.gz --output-type z --columns ID -a data/Homo_sapiens.r75.sort.vcf.gz data/MPRA_Variants.sort.vcf.gz
```

## Genomic sequences contexts from MPRA variants
Genomic sequences surrounding the assessed MPRA variants were extracted from the overlapping constructs and tested for any
mismatch or discontinuity in the overlapping sequences.
**Time consumming step, final files are provided**

```r
#Read constructs
constructs <- read.table(file = "data/TableS1_Constructs.tab",header = FALSE,sep = "\t")
names(constructs) <- c("variant","construct_pos","construct")

#Create data frame for final results and errors
df.final  <- NULL
df.errors <- NULL 

#Analyze each group of variants
for(variation in unique(constructs$variant)) {
  #Subset the 6 constructs per locus, i.e. 3 ref and 3 mut
  sub_const <- subset(constructs, constructs$variant == variation)
  
  #Create DNAStringSet objects from data.frames
  constr.fa  <- DNAStringSet(sub_const$construct)

  #Assign fasta headers to sequences
  names(constr.fa)  <- gsub(":" , "-", paste(sub_const$variant,sub_const$construct_pos, sep=";",collapse = NULL) )

  #Create names
  refname <- paste("Constructs", strsplit(names(constr.fa[1]),";")[[1]][1], "ref", sep=".",collapse = NULL)
  altname <- paste("Constructs", strsplit(names(constr.fa[4]),";")[[1]][1], "alt", sep=".",collapse = NULL)
  #Create fasta names
  refname.fa <- paste(paste("alignments/",refname, sep="",collapse = NULL), "fa", sep=".",collapse = NULL)
  altname.fa <- paste(paste("alignments/",altname, sep="",collapse = NULL), "fa", sep=".",collapse = NULL)
  #Create alignment names
  refname.aln <- paste(paste("alignments/",refname, sep="",collapse = NULL), "aln", sep=".",collapse = NULL)
  altname.aln <- paste(paste("alignments/",altname, sep="",collapse = NULL), "aln", sep=".",collapse = NULL)
  #Create pairwise alignment fasta name
  consensusname.fa <- paste("alignments/Consensus",strsplit(names(constr.fa[1]),";")[[1]][1], "fa",sep=".",collapse = NULL)
  #Create pairwise alignment fasta name
  consensusname.aln <- paste("alignments/Consensus",strsplit(names(constr.fa[1]),";")[[1]][1], "aln",sep=".",collapse = NULL)

    
  #Creata Fasta files for ref and alt
  writeXStringSet(constr.fa[c(1:3)],file= refname.fa, append = FALSE, format = "fasta", width = 145) 
  writeXStringSet(constr.fa[c(4:6)],file= altname.fa, append = FALSE, format = "fasta", width = 145) 
  
  #Perform mutiple alignments for each alt and ref Fasta file
  system( paste("clustalw2", refname.fa, sep=" ", collapse = NULL) )        
  system( paste("clustalw2", altname.fa, sep=" ", collapse = NULL) )        

  #Read aligments
  clustal.ref  <- readDNAMultipleAlignment( refname.aln, "clustal") 
  clustal.alt  <- readDNAMultipleAlignment( altname.aln, "clustal") 
  
  #Build consensus Matrix
  consMatrix.ref <- consensusMatrix(clustal.ref)
  consMatrix.alt <- consensusMatrix(clustal.alt)
  
  #Convert alignments to character arrays
  clustalChar.ref <- as.character(clustal.ref)
  clustalChar.alt <- as.character(clustal.alt)
  
  #Test sequence continuity for reference alignment
  valid <- check_seqcontinuity(clustalChar.ref)
  if(valid == 0){
    df.errors <- rbind(df.errors,variation)
    next
  }
  
  #Test sequence continuity for alternative alignment
  valid <- check_seqcontinuity(clustalChar.alt)
  if(valid == 0){
    df.errors <- rbind(df.errors,variation)
    next
  }

  #Get consensus sequence for ref and alt
  ref.consensus <- getConsensusSeq(consMatrix.ref)
  alt.consensus <- getConsensusSeq(consMatrix.alt)
  #Test return value
  if(class(ref.consensus) == "numeric" || class(alt.consensus) == "numeric"){
    df.errors <- rbind(df.errors,variation)
    next
  }
  
  #Arrange both consensus seq and identifiers in DF
  consensus.df <- data.frame(sequences = c(ref.consensus, alt.consensus),
                             identifiers = c(paste(strsplit(names(constr.fa[1]),";")[[1]][1], "ref", sep=".",collapse = NULL),
                                             paste(strsplit(names(constr.fa[4]),";")[[1]][1], "alt", sep=".",collapse = NULL)) )
  
  #Create DNAStringSet objects from data.frame
  consensus.fa  <- DNAStringSet(consensus.df$sequences)
  #Assign fasta headers to sequences
  names(consensus.fa) <- consensus.df$identifiers

  
  #Creata Fasta files for ref and alt consensus sequence
  writeXStringSet(consensus.fa, file= consensusname.fa, append = FALSE, format = "fasta", width = 145) 

  #Perform pairwise alignment for consensus sequences
  system( paste("clustalw2", consensusname.fa, sep=" ", collapse = NULL) )  
  
  
  #Extract variant information
  alleles <- strsplit(as.character(sub_const$variant[1]),":")[[1]][c(3,4)]
  alleles_size <- nchar(alleles)
  sizediff <- diff(alleles_size)
  
  #Read aligment
  clustal.con  <- readDNAMultipleAlignment(consensusname.aln, "clustal") 
  
  #Build consensus Matrix
  consMatrix.con <- consensusMatrix(clustal.con)
  
  ######################################################################
  ## Test if the assessed constructs are insertion, deletions or SNPs ##
  ## in order to test them more accurately.                           ##
  ######################################################################
  
  #SNP entry
  if(sizediff == 0){
    #Calculate indexes for splitting the consensus matrix in left 
    #sequence of variant, variant and right sequence of variant
    index_end_right_seq   <- dim(consMatrix.con)[2]
    index_end_left_seq    <- (index_end_right_seq - 1) /2
    index_start_var_seq   <- index_end_left_seq  + 1 
    index_end_var_seq     <- index_start_var_seq
    index_start_right_seq <- index_end_var_seq + 1
    
    #Subset consensus matrix
    consMatrix.left_seq  <- consMatrix.con[,1:index_end_left_seq]
    consMatrix.var_seq   <- consMatrix.con[,index_start_var_seq:index_end_var_seq]
    consMatrix.right_seq <- consMatrix.con[,index_start_right_seq:index_end_right_seq]
    
    #Convert alignment to character array
    clustalChar.con <- as.character(clustal.con)
    #Subset character array
    clustalChar.left_seq  <- substr(clustalChar.con,1,index_end_left_seq)
    clustalChar.var_seq   <- substr(clustalChar.con,index_start_var_seq,index_end_var_seq)
    clustalChar.right_seq <- substr(clustalChar.con,index_start_right_seq,index_end_right_seq)

    #Get consensus sequences for both flanked sequences
    left_seq.consensus  <- getConsensusSeq(consMatrix.left_seq)
    right_seq.consensus <- getConsensusSeq(consMatrix.right_seq)
    #Test return value
    if(class(left_seq.consensus) == "numeric" || class(right_seq.consensus) == "numeric"){
      df.errors <- rbind(df.errors,variation)
      next
    }
    
    #Ensure both sequences do not contain gaps
    if(grepl("-",clustalChar.con[1])  || grepl("-",clustalChar.con[2])  ){
      print(paste("A gap fas founded in SNP sequences at ",sub_const$variant,".",collapse = ""))
      df.errors <- rbind(df.errors,variation)
      next
    }
    #Ensure that both variant alleles are the same as those reported
    if(!(clustalChar.var_seq[1] == alleles[1] && clustalChar.var_seq[2] == alleles[2])){
      print(paste("Inconsistent alleles at ",sub_const$variant,".",collapse = ""))
      df.errors <- rbind(df.errors,variation)
      next
    }
    #Append processed entry to df.final, 1st alternative then reference allele
    var_info <- strsplit(as.character(sub_const$variant[1]),":")[[1]]
    df.final <- rbind(df.final,
                      c( var_info[1], as.numeric(var_info[2]) - 1, as.numeric(var_info[2]) -1 + alleles_size[1], "+", ".", 
                         "SNV", alleles[1], alleles[2], ".", 
                         paste0( c( tolower(left_seq.consensus),alleles[2],tolower(right_seq.consensus) ), collapse = "" ),
                         clustalChar.left_seq[2], clustalChar.right_seq[2],
                         nchar(left_seq.consensus), nchar(right_seq.consensus)))

    df.final <- rbind(df.final,
                      c( var_info[1], as.numeric(var_info[2]) - 1, as.numeric(var_info[2]) -1 + alleles_size[1], "+", ".", 
                         "SNV", alleles[1], alleles[1], ".", 
                         paste0( c( tolower(left_seq.consensus),alleles[1],tolower(right_seq.consensus) ), collapse = "" ),
                         clustalChar.left_seq[1], clustalChar.right_seq[1],
                         nchar(left_seq.consensus), nchar(right_seq.consensus)))
    
  #Insertion entry
  } else if( sizediff > 0){
    #Calculate indexes for splitting the consensus matrix in left 
    #sequence of variant, variant and right sequence of variant
    index_end_right_seq   <- dim(consMatrix.con)[2]
    index_end_left_seq    <- 96
    index_start_var_seq   <- index_end_left_seq  + 1 
    index_end_var_seq     <- index_start_var_seq + abs(sizediff)
    index_start_right_seq <- index_end_var_seq + 1
    
    #Subset consensus matrix
    consMatrix.left_seq  <- consMatrix.con[,1:index_end_left_seq]
    consMatrix.right_seq <- consMatrix.con[,index_start_right_seq:index_end_right_seq]
    
    #Convert alignment to character array
    clustalChar.con <- as.character(clustal.con)
    #Subset character array
    clustalChar.left_seq  <- substr(clustalChar.con,1,index_end_left_seq)
    clustalChar.var_seq   <- substr(clustalChar.con,index_start_var_seq,index_end_var_seq)
    clustalChar.right_seq <- substr(clustalChar.con,index_start_right_seq,index_end_right_seq)

    #Get consensus sequences for both flanked sequences
    left_seq.consensus  <- getConsensusSeq(consMatrix.left_seq)
    right_seq.consensus <- getConsensusSeq(consMatrix.right_seq)
    #Test return value
    if(class(left_seq.consensus) == "numeric" || class(right_seq.consensus) == "numeric"){
      df.errors <- rbind(df.errors,variation)
      next
    }

    #Ensure that both variant alleles are the same as those reported
    if(!( substr(clustalChar.var_seq[1],1,1) == alleles[1] && clustalChar.var_seq[2] == alleles[2] )){
      print(paste("Inconsistent alleles at ",sub_const$variant,".",collapse = ""))
      df.errors <- rbind(df.errors,variation)
      next
    }
    #Append processed entry to df.final, 1st alternative then reference allele
    var_info <- strsplit(as.character(sub_const$variant[1]),":")[[1]]
    df.final <- rbind(df.final,
                      c( var_info[1], as.numeric(var_info[2]) - 1, as.numeric(var_info[2]) - 1  + alleles_size[1], "+", ".", 
                         "insertion", alleles[1], alleles[2], ".", 
                         paste0( c( tolower(left_seq.consensus),alleles[2],tolower(right_seq.consensus) ), collapse = "" ),
                         gsub("-","",clustalChar.left_seq[2]), gsub("-","",clustalChar.right_seq[2]),
                        nchar(left_seq.consensus), nchar(right_seq.consensus)))

    df.final <- rbind(df.final,
                      c( var_info[1], as.numeric(var_info[2]) - 1, as.numeric(var_info[2]) - 1 + alleles_size[1], "+", ".", 
                         "insertion", alleles[1], alleles[1], ".", 
                         paste0( c( tolower(left_seq.consensus),alleles[1],tolower(right_seq.consensus) ), collapse = "" ),
                         gsub("-","",clustalChar.left_seq[1]), gsub("-","",clustalChar.right_seq[1]),
                        nchar(left_seq.consensus), nchar(right_seq.consensus)))
  #Deletion entry
  } else{
    #Calculate indexes for splitting the consensus matrix in left 
    #sequence of variant, variant and right sequence of variant
    index_end_right_seq   <- dim(consMatrix.con)[2]
    index_end_left_seq    <- 96
    index_start_var_seq   <- index_end_left_seq  + 1 
    index_end_var_seq     <- index_start_var_seq + abs(sizediff)
    index_start_right_seq <- index_end_var_seq + 1
    
    #Subset consensus matrix
    consMatrix.left_seq  <- consMatrix.con[,1:index_end_left_seq]
    consMatrix.right_seq <- consMatrix.con[,index_start_right_seq:index_end_right_seq]
    
    #Convert alignment to character array
    clustalChar.con <- as.character(clustal.con)
    #Subset character array
    clustalChar.left_seq  <- substr(clustalChar.con,1,index_end_left_seq)
    clustalChar.var_seq   <- substr(clustalChar.con,index_start_var_seq,index_end_var_seq)
    clustalChar.right_seq <- substr(clustalChar.con,index_start_right_seq,index_end_right_seq)

    #Get consensus sequences for both flanked sequences
    left_seq.consensus  <- getConsensusSeq(consMatrix.left_seq)
    right_seq.consensus <- getConsensusSeq(consMatrix.right_seq)
    #Test return value
    if(class(left_seq.consensus) == "numeric" || class(right_seq.consensus) == "numeric"){
      df.errors <- rbind(df.errors,variation)
      next
    }

    #Ensure that both variant alleles are the same as those reported
    if(!(clustalChar.var_seq[1] == alleles[1] && substr(clustalChar.var_seq[2],1,1) == alleles[2])){
      print(paste("Inconsistent alleles at ",sub_const$variant,".",collapse = ""))
      df.errors <- rbind(df.errors,variation)
      next
    }
    #Append processed entry to df.final, 1st alternative then reference allele
    var_info <- strsplit(as.character(sub_const$variant[1]),":")[[1]]
    df.final <- rbind(df.final,
                      c( var_info[1], as.numeric(var_info[2]) - 1, as.numeric(var_info[2]) - 1 + alleles_size[1], "+", ".", 
                         "deletion", alleles[1], alleles[2], ".", 
                         paste0( c( tolower(left_seq.consensus),alleles[2],tolower(right_seq.consensus) ), collapse = "" ),
                         gsub("-","",clustalChar.left_seq[2]), gsub("-","",clustalChar.right_seq[2]), 
                         nchar(left_seq.consensus), nchar(right_seq.consensus)))

    df.final <- rbind(df.final,
                      c( var_info[1], as.numeric(var_info[2]) - 1, as.numeric(var_info[2]) - 1 + alleles_size[1], "+", ".", 
                         "deletion", alleles[1], alleles[1], ".", 
                         paste0( c( tolower(left_seq.consensus),alleles[1],tolower(right_seq.consensus) ), collapse = "" ),
                         gsub("-","",clustalChar.left_seq[1]), gsub("-","",clustalChar.right_seq[1]),
                         nchar(left_seq.consensus), nchar(right_seq.consensus)))
    
  }
  
}

#Write  results
write.table(df.final[,c(-11,-12)], file = "data/MPRAs_reconstructed.mod.varSeq",sep = "\t",
            quote = FALSE ,row.names = FALSE, col.names = FALSE)
write.table(df.errors, file = "data/MPRAS.errors.tab",sep = "\t",
            quote = FALSE ,row.names = FALSE, col.names = FALSE)
```

##Variation-scan Analysis
### Sequences from MPRA variations and oligo-analysis background model


```bash


#convert-variations
convert-variations -from vcf -to varBed -i data/MPRA_Variants.annot75.vcf -o data/MPRA_Variants.annot75.varBed

#retrieve-variation-seq
retrieve-variation-seq -species Homo_sapiens -assembly GRCh37 -release 72 -mml 97 -format varBed \
-i data/MPRA_Variants.annot75.varBed -o data/MPRA_Variants.annot75.varSeq

#obtain sequences for oligo-analysis 
cut -f10 data/MPRA_Variants.annot75.varSeq | grep -v ";" | tail -n +2  > data/bg_sequences.txt

#oligo-analysis 
oligo-analysis -l 2 -i data/bg_sequences.txt -o data/MPRA_bgmodel.oa -format raw 

#Obtain all matrices for TFs relevant for erythrocytes
#Extract matrices from RSAT database
##Obtain names for matrices of interest
grep -i -E 'GATA1|KLF1|DHS|TAL1|ETS|FLI1|AP-1
|NFE2|LDB1' data/clusters_motif_names.tab | cut -f1 > data/matrixnames.list

##Extract matrices of interest from file
python3 source/matrices.py > data/mpra_matrices.txt

#variation-scan for 
variation-scan -i data/MPRA_Variants.annot75.varSeq -m data/mpra_matrices.txt -bg data/MPRA_bgmodel.oa  -o results/MPRA_casestudy.annot75.vs.nothresholds.table  -lth w_diff 0 -lth pval_ratio 1

#variation-scan for permuted matrices
variation-scan -i data/MPRA_Variants.annot75.varSeq -m data/mpra_permuted.tf -bg data/MPRA_bgmodel.oa  -o  results/MPRA_casestudy.annot75.vs.nothresholds.perm.table  -lth w_diff 0 -lth pval_ratio 1

```




## MPRA analysis

Taken from code supplied in Ulirsch et al., 2016.

This step is time consumoing, precomputed files are avialble in the folder *data/Precomputed*

### Read MPRA data
Read in raw barcode count data for MPRA:

```r
dir <- getwd() ####SET DIRECTORY but defult we use the local directoy where the Rmarkdown was started
MPRA <-
  data.frame(read_delim(
  file = paste0(dir, "/data/Raw/", "RBC_MPRA_minP_raw.txt"),
  delim = "\t",
  col_names = T,
  col_types = cols(chr = "c")
  ))
```


### Remove modified constructs (due to restriction enzyme cut sites):

```r
MPRA <- subset(MPRA, clean == "var")
```


### Combine plasmid counts from both replicates:

```r
MPRA <- as.data.frame(append(MPRA, list(K562_minP_DNA = MPRA$K562_minP_DNA1 + MPRA$K562_minP_DNA2),after = 13))
```


### Normalize counts to counts per million (CPM):

```r
for (i in 12:length(MPRA)) {
  MPRA[, i] <- (MPRA[, i] + 1) / 1000000 * sum(MPRA[, i])
  }
```

###Remove multi-allelic variants

```r
MPRA <- MPRA[!(MPRA$pos %in% c(2519416, 43842618, 46046134)), ]
#Exclude Duffy variant as well
MPRA <- MPRA[!(MPRA$pos %in% c(159174683)),]
```


###Log2(x+1) normalize and plot DNA barcode density

```r
MPRA[(length(MPRA) + 1):(length(MPRA) + length(MPRA[, 12:length(MPRA)]))] <- log(MPRA[, 12:(length(MPRA))], 2)
  
#Keep only barcodes with minimum log2(CPM+1) determined from density plot:
MPRA_minP <- MPRA[MPRA$K562_minP_DNA.1 >= 8,]
print(c("Percent of barcodes remaining: ", (dim(MPRA_minP) / dim(MPRA))[1]))
```

###Define activity [log2(mRNA/DNA) = log2(mRNA) - log2(DNA)]:

```r
attach(MPRA_minP)
MPRA_minP$K562_CTRL_RATIO_R1 <-
K562_CTRL_minP_RNA1.1 - K562_minP_DNA.1
MPRA_minP$K562_CTRL_RATIO_R2 <-
K562_CTRL_minP_RNA2.1 - K562_minP_DNA.1
MPRA_minP$K562_CTRL_RATIO_R3 <-
K562_CTRL_minP_RNA3.1 - K562_minP_DNA.1
MPRA_minP$K562_CTRL_RATIO_R4 <-
K562_CTRL_minP_RNA4.1 - K562_minP_DNA.1
MPRA_minP$K562_CTRL_RATIO_R5 <-
K562_CTRL_minP_RNA5.1 - K562_minP_DNA.1
MPRA_minP$K562_CTRL_RATIO_R6 <-
K562_CTRL_minP_RNA6.1 - K562_minP_DNA.1
MPRA_minP$K562_GATA1_RATIO_R1 <-
K562_GATA1_minP_RNA1.1 - K562_minP_DNA.1
MPRA_minP$K562_GATA1_RATIO_R2 <-
K562_GATA1_minP_RNA2.1 - K562_minP_DNA.1
MPRA_minP$K562_GATA1_RATIO_R3 <-
K562_GATA1_minP_RNA3.1 - K562_minP_DNA.1
MPRA_minP$K562_GATA1_RATIO_R4 <-
K562_GATA1_minP_RNA4.1 - K562_minP_DNA.1
detach(MPRA_minP)
```

###Quantile normalize activities, set median to 0, check, and combine replicates in separate DF:

```r
temp <- normalize.quantiles(as.matrix(MPRA_minP[, 38:47]))
temp <- temp - median(temp)
MPRA_minP.melt <-
melt(
MPRA_minP[c(
"chr",
"pos",
"ref",
"alt",
"type",
"bot",
"top",
"clean",
"oligo",
"construct",
"byallele",
"K562_CTRL_RATIO_R1",
"K562_CTRL_RATIO_R2",
"K562_CTRL_RATIO_R3",
"K562_CTRL_RATIO_R4",
"K562_CTRL_RATIO_R5",
"K562_CTRL_RATIO_R6",
"K562_GATA1_RATIO_R1",
"K562_GATA1_RATIO_R2",
"K562_GATA1_RATIO_R3",
"K562_GATA1_RATIO_R4"
)],
id = c(
"chr",
"pos",
"ref",
"alt",
"type",
"bot",
"top",
"clean",
"oligo",
"construct",
"byallele"
)
)
MPRA_minP.melt$value <- melt(temp[, 1:10])$value
```


###Derive activity estimates by collapsing barcodes and taking the median:

```r
CTRL.temp <- MPRA_minP.melt[grep("CTRL", MPRA_minP.melt$variable),]
CTRL.value <-
tapply(CTRL.temp$value, factor(CTRL.temp$byallele), median)
GATA1.temp <-
MPRA_minP.melt[grep("GATA1", MPRA_minP.melt$variable),]
GATA1.value <-
tapply(GATA1.temp$value, factor(GATA1.temp$byallele), median)
RATIO.temp <-
MPRA_minP.melt[!duplicated(MPRA_minP.melt$byallele),!(names(MPRA_minP.melt) %in% c("value", "variable"))]
MPRA_minP.ratio <-
merge(RATIO.temp,
data.frame(byallele = names(CTRL.value), CTRL.median = CTRL.value),
by = "byallele")
MPRA_minP.ratio <-
merge(MPRA_minP.ratio,
data.frame(byallele = names(GATA1.value), GATA1.median = GATA1.value),
by = "byallele")
```


###Create variable to indicate controls

```r
controls <-
  c("1 155271258",
  "X 55054634",
  "X 55054635",
  "X 55054636",
  "10 127505272")
  MPRA_minP.ratio <-
  as.data.frame(append(MPRA_minP.ratio, list(
  controls = ifelse(MPRA_minP.ratio$oligo %in% controls, 1, 0)
  ), after = 11))
  MPRA_minP.ratio <-
  MPRA_minP.ratio[!(MPRA_minP.ratio$construct %in% "1 159174683"), ]
```
  
###Determine active constructs and alleles with differences in activity

```r
  #Calculate p-values for overall enhancer activity
MPRA_minP.control.final <-
MPRA_minP.melt[grep("CTRL", MPRA_minP.melt$variable), ]
MPRA_minP.control.final.pvalues <- data.frame(matrix(ncol = 3))
colnames(MPRA_minP.control.final.pvalues) <-
c("byallele", "oligo", "Control_P")
for (j in levels(factor(MPRA_minP.control.final$byallele))) {
if (length(MPRA_minP.control.final[MPRA_minP.control.final$byallele == j, ]$value) >=
7) {
MPRA_minP.control.final.pvalues <-
rbind(
MPRA_minP.control.final.pvalues,
data.frame(
byallele = j,
oligo = MPRA_minP.control.final[MPRA_minP.control.final$byallele == j, ]$oligo[1],
Control_P = wilcox.test(
MPRA_minP.control.final[MPRA_minP.control.final$byallele == j, ]$value,
MPRA_minP.control.final[MPRA_minP.control.final$byallele != j, ]$value,
alternative = c("greater")
)$p.value
)
)
}
#if(j=="1 158459975 1/2 Ref") break
}

saveRDS(MPRA_minP.control.final.pvalues, file = "MPRA_minP.control.final.pvalues.rds")

#Determine alleles with differences in activity (2-sided Mann-Whitney-U test between alleles):

#Calculate p-values for differential activity between alleles
MPRA_minP.control.final <-
MPRA_minP.melt[grep("CTRL", MPRA_minP.melt$variable), ]
MPRA_minP.control.final.mut.pvalues <- data.frame(matrix(ncol = 3))
colnames(MPRA_minP.control.final.mut.pvalues) <-
c("construct", "oligo", "Control_P")
for (j in levels(factor(MPRA_minP.control.final$construct))) {
if (nlevels(factor(MPRA_minP.control.final[MPRA_minP.control.final$construct ==
j, ]$type)) == 2) {
MPRA_minP.control.final.mut.pvalues <-
rbind(
MPRA_minP.control.final.mut.pvalues,
data.frame(
construct = j,
oligo = MPRA_minP.control.final[MPRA_minP.control.final$byallele == j, ]$oligo[1],
Control_P = wilcox.test(value ~ factor(type), data = MPRA_minP.control.final[MPRA_minP.control.final$construct ==
j, ])$p.value
)
)
}
#if(j=="1 158459975 1/2 Ref") break
}

saveRDS(MPRA_minP.control.final.mut.pvalues, file = "MPRA_minP.control.final.mut.pvalues.rds")
```

###Read in already calculated p-values for activity and differential activity by allele:

```r
MPRA_minP.CTRL.pvalues <-
  readRDS(paste0(dir, "/data/Precomputed/","MPRA_minP.control.final.pvalues.rds"))

  MPRA_minP.GATA1.pvalues <-
  readRDS(paste0(dir, "/data/Precomputed/","MPRA_minP.GATA1.final.pvalues.rds"))

  MPRA_minP.CTRL.mut.pvalues <-
  readRDS(paste0(dir, "/data/Precomputed/","MPRA_minP.control.final.mut.pvalues.rds"))
  MPRA_minP.GATA1.mut.pvalues <-
  readRDS(paste0(dir, "/data/Precomputed/","MPRA_minP.GATA1.final.mut.pvalues.rds"))
```

###Merge test results and calculate FDR:

```r
  MPRA_minP.ratio <-
  merge(MPRA_minP.ratio,
  MPRA_minP.CTRL.pvalues,
  by = c("byallele", "oligo"))
  MPRA_minP.ratio <-
  merge(MPRA_minP.ratio,
  MPRA_minP.GATA1.pvalues,
  by = c("byallele", "oligo"))
  colnames(MPRA_minP.ratio)[15:16] <- c("CTRL.p", "GATA1.p")
  MPRA_minP.ratio$CTRL.padj <-
  qvalue(MPRA_minP.ratio$CTRL.p)[3]$qvalues
  MPRA_minP.ratio$GATA1.padj <-
  qvalue(MPRA_minP.ratio$GATA1.p)[3]$qvalues
  MPRA_minP.ratio <-
  merge(MPRA_minP.ratio,
  MPRA_minP.CTRL.mut.pvalues[, c(1, 3)],
  by = c("construct"))
  MPRA_minP.ratio <-
  merge(MPRA_minP.ratio,
  MPRA_minP.GATA1.mut.pvalues[, c(1, 3)],
  by = c("construct"))
  colnames(MPRA_minP.ratio)[19:20] <- c("CTRL.mut.p", "GATA1.mut.p")
  MPRA_minP.ratio$CTRL.mut.padj <-
  qvalue(MPRA_minP.ratio$CTRL.mut.p)[3]$qvalues
  MPRA_minP.ratio$GATA1.mut.padj <-
  qvalue(MPRA_minP.ratio$GATA1.mut.p)[3]$qvalues

#Read in tags
tags <-
  read.table(paste0(dir, "/data/Annotations/", "1000Genomes_Pilot3_LDpt8.tags"),
  sep = " ")
  tags$tagSNP <- do.call(paste, as.data.frame(cbind(tags$V1, tags$V2)))
  tags$oligo <- do.call(paste, as.data.frame(cbind(tags$V5, tags$V6)))
  colnames(tags) <-
  c(
  "tag_chr",
  "tag_pos",
  "tag_ref",
  "tag_alt",
  "chr",
  "pos",
  "ref",
  "alt",
  "tagSNP",
  "oligo"
  )
 #Calculate fold change
 REF <-
  MPRA_minP.ratio[MPRA_minP.ratio$type == "Ref", c("construct", "CTRL.median", "GATA1.median")]
  MUT <-
  MPRA_minP.ratio[MPRA_minP.ratio$type == "Mut", c("construct", "CTRL.median", "GATA1.median")]
  temp <- merge(MUT, REF, by = "construct")
  temp$CTRL.fc <- temp$CTRL.median.x - temp$CTRL.median.y
  temp$GATA1.fc <- temp$GATA1.median.x - temp$GATA1.median.y
  temp <- temp[, c("construct", "CTRL.fc", "GATA1.fc")]
  MPRA_minP.ratio <- merge(MPRA_minP.ratio, temp, by = "construct")
```

 

### MPRA functional variants and their effect sizes

```r
MPRA_minP.sig <-
  MPRA_minP.ratio[MPRA_minP.ratio$CTRL.mut.padj < 0.01 |
  MPRA_minP.ratio$GATA1.mut.padj < 0.01, ]
  MPRA_minP.sig.annot <- merge(MPRA_minP.sig, tags, by = "oligo")

  ###Write tables
write.table(MPRA_minP.ratio, file = paste0(dir, "/results/MPRA_results.tab"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(MPRA_minP.sig, file = paste0(dir, "/results/MPRA_sig.tab"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(MPRA_minP.sig.annot, file = paste0(dir, "/results/MPRA_sig.annot.tab"), sep="\t", quote=FALSE, row.names=FALSE)
```

##DeltaSVM

```bash

python3 source/VarSeq2Fasta.py data/MPRA_Variants.annot75.varSeq  data/MPRA_Variants.annot75
./deltasvm.pl data/MPRA_Variants.annot75.ref.fasta data/MPRA_Variants.annot75.alt.fasta data/SupplementaryTable_k562weights.txt results/MPRA_Variants_deltaSVM_output.txt

```




##Comparative analysis
###Read data 

```r
vs <-  read.table( paste0(dir, "/results/MPRA_casestudy.annot75.vs.nothresholds.table"), stringsAsFactors = FALSE, sep = "\t", header = TRUE, row.names = NULL, comment.char = ";") #Variation scan results with all matrices
tfnames <- read.table(paste0(dir, "/data/clusters_motif_names_no_duplicates.tab")) #Transcription factor name equivalency
colnames(tfnames) <- c("ac_motif", "tf_name")
mpra <- read.table(paste0(dir, "/results/MPRA_results.tab"), header = T, sep = "\t") #MPRA results

mpra_sig <- read.table(paste0(dir, "/results/MPRA_sig.annot.tab"), header = T, sep="\t") #Variants deemed as functional according to MPRA study

deltasvm_results <- read.table(paste0(dir, "/results/MPRA_Variants_deltaSVM_output.txt"))

deepsea_results <- read.csv(paste0(dir, "/results/MPRA_Variants_DeepSea.txt"))




## PERMUTED MATRICES VARIATION-SCAN RESULTS

vs_perm <- read.table(paste0(dir, "/results/MPRA_casestudy.annot75.vs.nothresholds.perm.table"), header = T, sep = "\t", comment.char = ";")



###Make "summary" of five permutations of each table 
coords <- str_match(as.character(vs_perm$var_coord), regex("((\\d|X)+).\\d+.(\\d+)"))
vs_perm$posconcat <- paste(coords[,2], coords[,4], sep = ":")
vs_perm$clustervar <- paste(vs_perm$posconcat, substr(vs_perm$ac_motif, 1, str_length(vs_perm$ac_motif)-6))
pvalratio_sum <- (tapply(vs_perm$pval_ratio, vs_perm$clustervar, mean))
pvaldf <- data.frame(clustervar = names(pvalratio_sum), pval_ratio = pvalratio_sum)
bestw_sum <- tapply(vs_perm$best_w, vs_perm$clustervar, mean)
bestwdf <- data.frame(clustervar = names(bestw_sum), best_w = bestw_sum)
worstw_sum <-tapply(vs_perm$worst_w, vs_perm$clustervar, mean)
worstwdf <- data.frame(clustervar = names(worstw_sum), worst_w = worstw_sum)

vs_permuted_sum <- merge(pvaldf, bestwdf)
vs_permuted_sum <- merge(vs_permuted_sum, worstwdf)
temp <- data.frame(posconcat = vs_perm$posconcat, clustervar = vs_perm$clustervar)
temp <- unique(temp)
vs_permuted <- merge(temp, vs_permuted_sum)

##ANNOTATE TABLES

#Define "posconcat" variable in dataframes for all results, in order to be able to compare variants between both studies. This variable concatenates the chromosome and position of each variant

coords <- str_match(as.character(vs$var_coord), regex("((\\d|X)+).\\d+.(\\d+)"))
vs$posconcat <- paste(coords[,2], coords[,4], sep = ":")

mpra_sig$posconcat <- paste(mpra_sig$chr.x, mpra_sig$pos.x, sep = ":")

mpra$posconcat <- paste(mpra$chr, mpra$pos, sep = ":")

mpra_nonsig <- mpra[!mpra$oligo %in% mpra_sig$oligo,] #Variables in MPRA not deemed as functional 
```



Filter to extract results with a sign change in the weight, and a fold change in p-value.

```r
vs_sign_change <- vs[vs[,"best_w"]/vs[,"worst_w"]<=0,]
vs_sign_change <- merge(vs_sign_change, tfnames)
vs_fold_change <- vs_sign_change[vs_sign_change[,"pval_ratio"]>=10,]
vs_fold_change_pval <- vs_fold_change[vs_fold_change[,"best_pval"]<=0.0001,]

vs_perm_sign_change <- vs_perm[vs_perm[,"best_w"]/vs_perm[,"worst_w"]<=0,]
vs_perm_fold_change <- vs_perm_sign_change[vs_perm_sign_change[,"pval_ratio"]>=10,]
vs_perm_fold_change_pval <- vs_perm_fold_change[vs_perm_fold_change[,"best_pval"]<=0.0001,]
```


##UpsetR figure

```r
vs_change <- vs_fold_change_pval$posconcat
mpra_change <- unique(mpra_sig$posconcat)
vs_nochange <- setdiff(mpra$posconcat, unique(vs_fold_change_pval$posconcat))
mpra_nochange <- unique(mpra_nonsig$posconcat)
rslist <- list(variation_scan_change = vs_change, 
               MPRA_change =  mpra_change, 
               variation_scan_nochange =  vs_nochange, 
               MPRA_nochange =  mpra_nochange)

pdf("UpsetR_exclude_corrected.pdf")
upset(fromList(rslist), order.by = "freq" )
dev.off()
```

###UpsetR figure for permuted matrices

```r
vs_change <- vs_perm_fold_change_pval$posconcat
mpra_change <- unique(mpra_sig$posconcat)
vs_nochange <- setdiff(mpra$posconcat, unique(vs_perm_fold_change_pval$posconcat))
mpra_nochange <- unique(mpra_nonsig$posconcat)
rslist <- list(variation_scan_change = vs_change, 
               MPRA_change =  mpra_change, 
               variation_scan_nochange =  vs_nochange, 
               MPRA_nochange =  mpra_nochange)


pdf("UpsetR_exclude_perm_corrected.pdf")
upset(fromList(rslist), order.by = "freq" )
dev.off()
```


###ROC Curve for permuted matrices and real matrices

```r
vs_perm_pvalratio_summary <- tapply(vs_perm$pval_ratio, vs_perm$clustervar, mean)
vs_perm_clustervarmean <- tapply(vs_perm$pval_ratio, vs_perm$clustervar, mean)
vs_perm_clustervarmean_df <- data.frame(clustervar = names(vs_perm_clustervarmean), clustervarmean = vs_perm_clustervarmean)
vs_perm_subset <- data.frame(vs_perm$posconcat,  vs_perm$pval_ratio, vs_perm$clustervar)
colnames(vs_perm_subset) <- c("posconcat", "pval_ratio", "clustervar")
vs_perm_subset <- merge(vs_perm_subset, vs_perm_clustervarmean_df)

vs_perm_pval <- tapply(vs_perm_subset$clustervarmean, vs_perm_subset$posconcat, max)
vs_pval <- tapply(vs$pval_ratio, vs$posconcat, max)


vs_df <- data.frame(posconcat = names(vs_pval), pval = vs_pval)
vs_df$ispositive <- as.numeric(vs_df$posconcat %in% mpra_sig$posconcat)
vs_df <-  vs_df[order(vs_df$pval, decreasing = T),]

vs_perm_df <- data.frame(posconcat = names(vs_perm_pval), pval = vs_perm_pval)
vs_perm_df$ispositive <- as.numeric(vs_perm_df$posconcat %in% mpra_sig$posconcat)

separate <- strsplit(as.character(deltasvm_results$V1), ":")
chromosome <- unlist(separate)[2*(1:length(separate))-1]
pos <- unlist(separate)[2*(1:length(separate))]
deltasvm_results$posconcat <- paste(chromosome, (as.numeric(pos)+1), sep = ":")
deltasvm_results$ispositive <- as.numeric(deltasvm_results$posconcat %in% mpra_sig$posconcat)
deltasvm_results <- deltasvm_results[order(deltasvm_results$V2, decreasing = T),]
vs_perm_df <- vs_perm_df[order(vs_perm_df$pval, decreasing = T),]


deepsea_results$posconcat <- paste(substring(deepsea_results$chr,4), deepsea_results$pos, sep = ":")
deepsea_results$ispositive <- as.numeric(deepsea_results$posconcat %in% mpra_sig$posconcat)
deepsea_results <- deepsea_results[order(deepsea_results$Functional.significance.score, decreasing = T),]


roc_obj1 <- roc(vs_perm_df$ispositive, vs_perm_df$pval, plot = T, print.auc = T)
roc_obj2 <- roc(deltasvm_results$ispositive, deltasvm_results$V2, plot = T, print.auc = T)
roc_obj3 <- roc(vs_df$ispositive, vs_df$pval, plot = T, print.auc = T)
roc_obj4 <- roc(deepsea_results$ispositive, deepsea_results$Functional.significance.score, plot = T, print.auc = T)

#Calculate specificity at threshold p-values

tp <- vs_df$ispositive
fn <- rev(tp)
limit_10 <- which.min(abs(vs_df$pval-10)^2)
limit_100 <- which.min(abs(vs_df$pval-100)^2)
limit_1000 <- which.min(abs(vs_df$pval-1000)^2)



recall_10 <- if (sum(tp[1:limit_10]) + sum(fn[limit_10:length(fn)]) == 0) 0 else (sum(tp[1:limit_10]) / (sum(tp[1:limit_10]) + sum(fn[limit_10:length(fn)])))
recall_100 <- if (sum(tp[1:limit_100]) + sum(fn[limit_100:length(fn)]) == 0) 0 else (sum(tp[1:limit_100]) / (sum(tp[1:limit_100]) + sum(fn[limit_100:length(fn)])))
recall_1000 <- if (sum(tp[1:limit_1000]) + sum(fn[limit_1000:length(fn)]) == 0) 0 else (sum(tp[1:limit_1000]) / (sum(tp[1:limit_1000]) + sum(fn[limit_1000:length(fn)])))






pdf("Tool_Comparison_ROC.pdf")
ggroc(list(permuted = roc_obj1, deltasvm = roc_obj2, variationscan = roc_obj3, deepsea =roc_obj4))+
  coord_equal() +
  theme(aspect.ratio = 1) + 
  geom_hline(mapping = aes(x = 0, xend = 1, y = 0, yend = 1), yintercept = recall_10, color = "darkgrey", linetype = "dashed") +
  geom_hline(mapping = aes(x = 0, xend = 1, y = 0, yend = 1), yintercept = recall_100, color = "darkgrey", linetype = "dashed") +
  geom_hline(mapping = aes(x = 0, xend = 1, y = 0, yend = 1), yintercept = recall_1000, color = "darkgrey", linetype = "dashed") 
dev.off()
```
