
## 1) Obtain ChIPseq peaks from literature

In [1] the authors describe the genomic targets of transcription factor VRN1 in barley (Hordeum vulgare),
which they map to individual contigs of the 2012 reference genome (v1) [2].

The following supplementary files were downloaded, which contain peaks and neighbor genes:

+ ncomms6882-s3.xlsx -> converted to ncomms6882-s3.bed with cut -f 2,3,4,8,9,12 

+ ncomms6882-s4.xlsx -> DE genes in RNAseq experiment (near peaks)

## 2) check peaks contain the Vrn1 motif and can be matched in 2017 genome (v2) 

The following operations require the v3 barley genome [3] and BLAST+ [4]:

```
perl _check_peaks_coords.pl > ncomms6882-s3.2017.bed

sort -k1,1 -k2,2n -k3,3n -k4,4n  ncomms6882-s3.2017.bed > ncomms6882-s3.2017.sort.bed
```

## 3) intersect peaks and barley variants from Ensembl Plants (EG42)

This step uses bedtools [5] and Ensembl Plants [6]:

```
bedtools --version 
#bedtools v2.26.0

bedtools intersect -sorted -a hordeum_vulgare.vcf -b ncomms6882-s3.2017.sort.bed -wa > hordeum_vulgare_ncomms6882-s3.2017.vcf

wc hordeum_vulgare_ncomms6882-s3.2017.vcf
#1604
```

## 4) Analyze overlapping variants in RSAT

The following steps require RSAT, either a local installation or the web server at http://plants.rsat.eu, 
and the Vrn1 motif obtained from <http://floresta.eead.csic.es/footprintdb/index.php?motif=AY750993:VRN1:EEADannot> [7]:

```
$RSAT/perl-scripts/convert-variations -v 0  -i hordeum_vulgare_ncomms6882-s3.2017.vcf -from vcf -to varBed

$RSAT/bin/retrieve-variation-seq -source plants -species Hordeum_vulgare -assembly IBSCv2 \
	-i hordeum_vulgare_ncomms6882-s3.2017.varBed -format varBed -mml 30

$RSAT/bin/variation-scan -v 1 -m VRN1.tf -m_format transfac -i hordeum_vulgare_ncomms6882-s3.2017.varSeq \
	-bg 2nt_upstream-noorf_Hordeum_vulgare.IBSCv2.42-ovlp-1str.freq.gz \
	-lth score 1  -lth w_diff 1  -lth pval_ratio 10  -uth pval 1e-3 
```

## 5) Retrieve genes linked to these variants

```
cut -f 2 hordeum_vulgare_ncomms6882-s3.2017.variants-seq_result | grep vcZ > \
	hordeum_vulgare_ncomms6882-s3.2017.variant.list

cat hordeum_vulgare_ncomms6882-s3.2017.variant.list | sort -u | wc
#13

grep -P '^#' hordeum_vulgare.vcf > hordeum_vulgare_ncomms6882-s3.2017.variant.vcf
grep -f hordeum_vulgare_ncomms6882-s3.2017.variant.list hordeum_vulgare.vcf >> \
	hordeum_vulgare_ncomms6882-s3.2017.variant.vcf

bedtools intersect -sorted -a ncomms6882-s3.2017.sort.bed -wa \
	-b hordeum_vulgare_ncomms6882-s3.2017.variant.vcf > hordeum_vulgare_ncomms6882-s3.2017.variant.annotation

grep -c MLOC ncomms6882-s4.tsv
#38

cut -f 4 hordeum_vulgare_ncomms6882-s3.2017.variant.annotation > MLOC.list
	
grep -f MLOC.list ncomms6882-s4.tsv > hordeum_vulgare_ncomms6882-s3.2017.variant.MLOC

cat hordeum_vulgare_ncomms6882-s3.2017.variant.MLOC
#AK371349	117	53	2,19	Peak_284	31	morex_contig_63135	MLOC_73196	Amino acid permease
#MLOC_79452	76	135	0,56	Peak_036	112	morex_contig_85083	MLOC_79452	Unknown protein
```

This means that 2 genes (out of 38) for which the expression is known to be affected by VRN1 binding harbour natural variants
in VRN1 sites.


## Literature

[1] Deng et al (2015) Direct links between the vernalization response and other key traits of cereal crops. Nature Communications volume 6, Article number: 5882 (2015) <https://www.nature.com/articles/ncomms6882>

[2] IBSC (2012) A physical, genetic and functional sequence assembly of the barley genome. Nature. 491:711-16. doi:10.1038/nature11543

[3] Mascher et al. 2017 A chromosome conformation capture ordered sequence of the barley genome. Nature. 544:427-433. doi:10.1038/nature2204

[4] Camacho et al. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.

[5] Quinlan & Hall (2010) BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26(6):841–842 <https://doi.org/10.1093/bioinformatics/btq033>

[6] Kersey et al (2018) Ensembl Genomes 2018: an integrated omics infrastructure for non-vertebrate species. Nucleic Acids Research 46(D1):D802–D808 <https://doi.org/10.1093/nar/gkx1011>

[7] Sebastian & Contreras-Moreira (2014) footprintDB: a database of transcription factors with annotated cis elements and binding interfaces. Bioinformatics 30, 258–65 <https://academic.oup.com/bioinformatics/article/30/2/258/226399>
