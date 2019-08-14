use strict;
use warnings;

## data files

# converted from https://media.nature.com/original/nature-assets/ncomms/2015/150106/ncomms6882/extref/ncomms6882-s3.xlsx,
# see README.md
my $bedfile = 'ncomms6882-s3.bed';

# from ftp://ftpmips.helmholtz-muenchen.de/plants/barley/public_data/sequences/assembly3_WGSMorex_renamed_blastable_carma.zip
my $oldContigsFasta = 'assembly3_WGSMorex_renamed_blastable_carma.fasta';

# from ftp://ftp.ensemblgenomes.org/pub/release-43/plants/fasta/hordeum_vulgare/dna/Hordeum_vulgare.IBSC_v2.dna_rm.toplevel.fa.gz
my $newGenomeFasta = 'Hordeum_vulgare.IBSC_v2.dna_rm.toplevel.fa';

## binaries, not supplied, should be installed in your system

my $bedtoolsEXE = 'bedtools getfasta';
my $blastnEXE = 'ncbi-blast-2.9.0+/bin/blastn';

##########################################

my ($chr,$start,$end,$mloc,$strand);
my ($match,$mchr,$mstart,$mend); # matched coords
my ($seq,$new_seq,$site,$rcsite);

open(BED,"<",$bedfile) ||
	die "# cannot read $bedfile\n";
while(<BED>){
	#morex_contig_157888	76	1100	AAAGAAATGG		
   #morex_contig_2547616	4362	5554	ACCAAAAAAG	MLOC_36994	-
	next if(/^#/);
	chomp;
	my @data = split(/\t/,$_);

   $chr = $data[0];
	$start = $data[1];
	$end = $data[2];
	$site = $data[3];
   $mloc = $data[4] || 'NA';
   
   # save BED format
	open(TMPBED,">_tmp.bed");
	print TMPBED "$chr\t$start\t$end\n";
	close(TMPBED);

	$rcsite = $site;
	$rcsite =~ tr/acgtACGT/tgcaTGCA/;
   $rcsite = reverse($rcsite);

   # get original peak sequence
   $seq = `$bedtoolsEXE -fi $oldContigsFasta -bed _tmp.bed -tab`;
	$seq = (split(/\s+/,$seq))[1];

   # check fragment contains VRN1 site
	$match = 0;
	if($seq =~ m/$site/ || $seq =~ m/$rcsite/){ $match = 1 }

   # save FASTA format
   open(TMPFA,">_tmp.fa");
   print TMPFA ">$chr $start $end\n$seq\n";
   close(TMPFA);

	# find perfect match in new (2017) genome
	my @blast;
	$new_seq = '';
	open(BLASTN,"$blastnEXE -db $newGenomeFasta -outfmt \"6 saccver pident mismatch gapopen qstart qend sstart send evalue bitscore sseq\" -query _tmp.fa |");
	while(<BLASTN>){
		#chr4H	100.000	0	0	227	441	147267501	147267715	3.37e-108	398	AGGTATA...
		chomp;
		@blast = split(/\t/,$_);
		next if($blast[1] < 100 || $blast[2] > 0 || $blast[3] > 0);
		$new_seq = $blast[10];
		last;
	}
	close(BLASTN);

   next if($new_seq eq '');

   # check matched (new) fragment contains VRN1 site
   if($new_seq =~ m/$site/ || $new_seq =~ m/$rcsite/){ $match++ }

   next if($match < 2);

	# get BLASTN match coords and reverse comp strands
   ($mchr,$mstart,$mend) = ($blast[0],$blast[6],$blast[7]);
	if($mstart>$mend){ ($mstart,$mend) = ($blast[7],$blast[6]) }

	printf("%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%d\t%s\t%s\n",
		$mchr,$mstart,$mend,$mloc,
		$chr,$start,$end,
		$site,$rcsite,$match,
		$new_seq,$seq);
}
close(BED);
