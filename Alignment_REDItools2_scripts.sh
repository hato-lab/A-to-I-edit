
############# Genome Setup

mkdir -p ./references/STAR/GRCh38
cd references/STAR/GRCh38/
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.primary_assembly.annotation.gtf.gz
gunzip *
module load samtools/1.9
samtools faidx GRCh38.primary_assembly.genome.fa
rm chrfile.txt; for i in $(seq 1 22) X Y M; do echo chr${i} >>chrfile.txt; done
remove_ids=($(cat chrfile.txt));echo ${remove_ids[@]}
samtools faidx -o GRCh38.primary_assembly.genome.filtered.fa GRCh38.primary_assembly.genome.fa "${remove_ids[@]}"

sed -i '0,/chr1/s//chr1 1/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr2/s//chr2 2/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr3/s//chr3 3/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr4/s//chr4 4/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr5/s//chr5 5/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr6/s//chr6 6/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr7/s//chr7 7/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr8/s//chr8 8/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr9/s//chr9 9/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr10/s//chr10 10/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr11/s//chr11 11/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr12/s//chr12 12/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr13/s//chr13 13/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr14/s//chr14 14/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr15/s//chr15 15/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr16/s//chr16 16/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr17/s//chr17 17/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr18/s//chr18 18/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr19/s//chr19 19/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr20/s//chr20 20/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr21/s//chr21 21/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chr22/s//chr22 22/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chrX/s//chrX X/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chrY/s//chrY Y/' GRCh38.primary_assembly.genome.filtered.fa
sed -i '0,/chrM/s//chrM MT/' GRCh38.primary_assembly.genome.filtered.fa

###### genomeGenerate

STAR --runMode genomeGenerate \
--genomeDir AtoIPipeline/Human/references/STAR/GRCh38/ \
--genomeFastaFiles AtoIPipeline/Human/references/STAR/GRCh38/GRCh38.primary_assembly.genome.filtered.fa \
--sjdbGTFfile AtoIPipeline/Human/references/STAR/GRCh38/gencode.v41.primary_assembly.annotation.gtf \
--runThreadN 16

###### substituted genome

mkdir -p AtoIPipeline/Human/references/STAR/GRCh38.a2g/ 
        #GRCh38.primary_assembly.genome.{substitution}.fa;
        <AtoIPipeline/Human/references/STAR/GRCh38/GRCh38.primary_assembly.genome.filtered.fa AtoIPipeline/Human/references/xsubstituteFASTX.sh a2g FASTA > AtoIPipeline/Human/references/STAR/GRCh38.a2g/GRCh38.primary_assembly.genome.a2g.fa
        # A fake chromosome containing some of the missing nucleotides seems to be necessary for STAR indexing to complete in a sane amount of time.
        
        printf ">Fake\n" >> AtoIPipeline/Human/references/STAR/GRCh38.a2g/GRCh38.primary_assembly.genome.a2g.fa \
        for one in {A,C,T,G}; do for two in {A,C,T,G}; do for three in {A,C,T,G}; do for four in {A,C,T,G}; do printf "${one}${two}${three}${four}" >> AtoIPipeline/Human/references/STAR/GRCh38.a2g/GRCh38.primary_assembly.genome.a2g.fa; done; done; done; done
        printf "\n" >> AtoIPipeline/Human/references/STAR/GRCh38.a2g/GRCh38.primary_assembly.genome.a2g.fa
        
###### genomeGenerate substituted

STAR --runMode genomeGenerate \
--genomeDir AtoIPipeline/Human/references/STAR/GRCh38.a2g/ \
--genomeFastaFiles AtoIPipeline/Human/references/STAR/GRCh38.a2g/GRCh38.primary_assembly.genome.a2g.filtered.fa \
--sjdbGTFfile AtoIPipeline/Human/references/STAR/GRCh38/gencode.v41.primary_assembly.annotation.gtf \
--runThreadN 16

###### First alignment

 # Run STAR in alignment mode. # zcat is necessary when input read files are gzipped. # RAM limit should be in bytes.
 #Running for each sample (names in config.yaml)
mkdir -p output/firstSTARalignment/; 
array=($(tail -n 20 config.yaml))
for sample in $(seq 0 19); 
do str1=($(ls FASTQ/${array[$sample]}*1.fq.gz)); str2=($(ls FASTQ/${array[$sample]}*2.fq.gz)); \
       
        STAR --runMode alignReads \
        --runThreadN 12 \
        --runDirPerm All_RWX \
        --genomeDir ./references/STAR/GRCh38/ \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesType Fastx \
        --readFilesIn $str1 $str2 \
        --readFilesCommand zcat \
        --limitBAMsortRAM 600000000000 \
        --outFileNamePrefix output/firstSTARalignment/${array[$sample]}. \
        --outReadsUnmapped Fastx \
        --outMultimapperOrder Random \
        --outBAMsortingThreadN 12 \
done

####### convert_unmapped_read_bases
 
mkdir -p output/substitutedUnmappedReads/ 

for sample in $(seq 0 19); do for pair in 1 2; do
<output/firstSTARalignment/${array[$sample]}.Unmapped.out.mate${pair} ./xsubstituteFASTX.sh a2g FASTQ | gzip > output/substitutedUnmappedReads/${array[$sample]}.Unmapped.out.mate${pair}.a2g.fastq.gz; done; done

#align unmapped reads

mkdir -p output/secondSTARalignment/
array=($(tail -n 20 config.yaml))
for sample in $(seq 0 19);
do str1=output/substitutedUnmappedReads/${array[$sample]}.Unmapped.out.mate1.a2g.fastq.gz; str2=output/substitutedUnmappedReads/${array[$sample]}.Unmapped.out.mate2.a2g.fastq.gz;   
        STAR --runMode alignReads \
        --runThreadN 12 \
        --runDirPerm All_RWX \
        --genomeDir references/STAR/GRCh38.a2g \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesType Fastx \
        --readFilesIn $str1 $str2 \
        --readFilesCommand zcat \
        --limitBAMsortRAM 600000000000 \
        --outFileNamePrefix output/secondSTARalignment/${array[$sample]}. \
        --outReadsUnmapped None \
        --outMultimapperOrder Random \
        --outBAMsortingThreadN 12 ;
done

###### xgetOrigRead.awk

#!/usr/bin/env gawk

# This script should take three arguments.
# The first two arguments should be fastq files corresponding to the
# first and second mate read pairs. They should have exactly 4 lines per
# read. The third argument should be a SAM file.

# Define a function for reverse complementing.
function revcomp(arg) {
	o = ""
	for(i = length(arg); i > 0; i--)
		o = o c[substr(arg, i, 1)]
	return(o)
}
# Establish the dictionary array for reverse complementing.
BEGIN {
	c["A"] = "T"; 
	c["C"] = "G"; 
	c["G"] = "C"; 
	c["T"] = "A";
	c["a"] = "t"; 
	c["c"] = "g"; 
	c["g"] = "c"; 
	c["t"] = "a";
	c["N"] = "N";
	c["n"] = "n"
}
# While reading in the first fastq file, make an array.
FILENAME==ARGV[1]  {
	# Grab the read name, excluding the starting @ symbol.
	readname = substr($1, 2)
	# Use getline to pull the next line into a variable.
	getline readseq
	# Add the read to an array based on the name.
	firstmate[readname] = readseq
	# Jump to the 4th line of the FASTQ read.
	getline
	getline
	# Ignore the rest of the patterns for this line.
	next
}
# Do this all again for the second read mates.
FILENAME==ARGV[2]  {
	# Grab the read name, excluding the starting @ symbol.
	readname = substr($1, 2)
	# Use getline to pull the next line into a variable.
	getline readseq
	# Add the read to an array based on the name.
	secondmate[readname] = readseq
	# Jump to the 4th line of the FASTQ read.
	getline
	getline
	# Ignore the rest of the patterns for this line.
	next
}
# Now read in the SAM file and put out the correct read sequence.
# You can check https://samtools.github.io/hts-specs/SAMv1.pdf section 1.4 for info on SAM columns and flags.
# Print the header lines unchanged.
$0 ~ /^@/ {print}
# Correct the read sequences.
$0 !~ /^@/ {
	readname = $1
	# Check if first or second read using the flag column $2.
	# Use bitwise rightshift to check flag 0x40
	if (and(rshift($2, 6), 1) == 1) {
		$10 = firstmate[readname]
	} else {
		$10 = secondmate[readname]
	}
	# Check if the read is mapped to the minus strand with flag 0x10.
	if (and(rshift($2, 4), 1) == 1) {
		$10 = revcomp($10)
	} 
	# Output the modified line.
	print $0
}


###### unsubstitute_bam_files
   
        module load samtools/1.10; 
        mkdir -p output/unsubstituted_second_alignments/
        array=($(tail -n 20 config.yaml))
        for sample in $(seq 0 19); 
        # Use a GNU-awk script to read in the original read files and the substituted BAM files, and output the correct read sequences.
        # Use samtools v1.10+ and process substitution to uncompress the BAM file on the fly.
        do first_antisense_mates=output/firstSTARalignment/${array[$sample]}.Unmapped.out.mate1;
        second_sense_mates=output/firstSTARalignment/${array[$sample]}.Unmapped.out.mate2;
        substituted_bam=output/secondSTARalignment/${array[$sample]}.Aligned.sortedByCoord.out.bam;
        output=output/unsubstituted_second_alignments/${array[$sample]}.a2g.Unsubstituted.Aligned.sortedByCoord.out.bam;
        gawk -v OFS='\t' -f xgetOrigRead.awk ${first_antisense_mates} ${second_sense_mates} <(samtools view -h ${substituted_bam}) | samtools view -b > ${output} ;
        done
        
####### merge_bam_files
        
        module load samtools/1.10;
        mkdir -p output/merged_alignments/;
        # Merge using samtools v1.10+.
        array=($(tail -n 20 config.yaml))
        for sample in $(seq 0 19); 
        do first_bam=output/firstSTARalignment/${array[$sample]}.Aligned.sortedByCoord.out.bam;
        unsubstituted_second_bam=output/unsubstituted_second_alignments/${array[$sample]}.a2g.Unsubstituted.Aligned.sortedByCoord.out.bam;
        merged_bam=output/merged_alignments/${array[$sample]}.a2g.Merged.Aligned.sortedByCoord.out.bam;
        samtools merge --threads 12 ${merged_bam} ${first_bam} ${unsubstituted_second_bam};
        done

   
######index_fasta:

   module load samtools/1.10; 
   samtools faidx references/STAR/GRCh38/GRCh38.primary_assembly.genome.fa;
   samtools faidx references/STAR/GRCh38.a2g/GRCh38.primary_assembly.genome.a2g.fa;    
   
######index_bam:
    
   module load samtools/1.10; 
   array=($(tail -n 20 config.yaml))
   for sample in $(seq 0 19); 
   do merged_bam=output/merged_alignments/${array[$sample]}.a2g.Merged.Aligned.sortedByCoord.out.bam;
   samtools index -@ 12 $merged_bam;
   done

#######											#######
## From https://github.com/BioinfoUNIBA/REDItools2   ##
#######											#######

######get_coverage_file
		
        for sample in $(seq 0 19); do array=($(tail -n 20 config.yaml));
			mkdir -p output/coverage/${array[$sample]}.a2g;
			REDItools2/extract_coverage.sh output/merged_alignments/${array[$sample]}.a2g.Merged.Aligned.sortedByCoord.out.bam output/coverage/${array[$sample]}.a2g/ references/STAR/GRCh38/GRCh38.primary_assembly.genome.fa.fai; done

######run_REDItools2

	for sample in $(seq 0 19); do export LD_LIBRARY_PATH="/N/soft/rhel7/python/2.7.16/lib/libfabric:$LD_LIBRARY_PATH"; 
	mpirun -n 12 python REDItools2/src/cineca/parallel_reditools.py -f output/merged_alignments/${array[$sample]}.a2g.Merged.Aligned.sortedByCoord.out.bam \
	-o output/reditools2/${array[$sample]}.a2g.reditools_table.txt.gz -r references/STAR/GRCh38/GRCh38.primary_assembly.genome.fa \
	-t output/temp/${array[$sample]}.a2g/ -Z references/STAR/GRCh38/GRCh38.primary_assembly.genome.fa.fai -G output/coverage/${array[$sample]}.a2g/${array[$sample]}.a2g.Merged.Aligned.sortedByCoord.out.cov \
	-D output/coverage/${array[$sample]}.a2g/ -q 30 -s 2 -C; done

######merge_REDItools2

bash xREDrun xMergeRedItools 
for sample in $(seq 0 19); do REDItools2/merge.sh output/temp/${array[$sample]}.a2g output/reditools2/${array[$sample]}.a2g.reditools_table.txt.gz 6; done

######find_every_edited_site

mkdir completeOutput/a2g_sites
array=($(cat completeOutput/samples.txt))
for file in $(seq 0 59); 
do echo ${array[$file]}; zcat completeOutput/reditools2/${array[$file]}.a2g.reditools_table.txt.gz| awk -v FS="\\t" -v OFS="\\t" ' $8 == "AG" {{ print $1, $2 - 1, $2, $1 ":" $2 "-" $2 }} {{fflush()}} ' | sort --buffer-size=12000M --parallel=1 -u -k1.1,1.1 -k1,1V -k2,2g  > completeOutput/a2g_sites/${array[$file]}_a2g_site.bed;
done

for file in $(seq 0 59); 
do zcat output/reditools2/${array[$file]}.a2g.reditools_table.txt.gz| awk -v FS="\\t" -v OFS="\\t" ' $8 == "AG" {{ print $1, $2 - 1, $2, $1 ":" $2 "-" $2 }} {{fflush()}} ' | sort --buffer-size=12000M --parallel=1 -u -k1.1,1.1 -k1,1V -k2,2g  > output/a2g_sites/${array[$file]}_a2g_site.bed; done

cat completeOutput/a2g_sites/* > completeOutput/a2g_sites/every_a2g_site.bed


######count_each_substitution

module load tabix/1.10;
printf "Region\tPosition\tReference\tStrand\tCoverage-q30\tMeanQ\tBaseCount[A,C,G,T]\tAllSubs\tFrequency\tgCoverage-q30\tgMeanQ\tgBaseCount[A,C,G,T]\tgAllSubsgFrequency\tSample\n" > completeOutput/a2g_sites/all_a2g.reditools_table.txt; 
for file in $( echo completeOutput/reditools2/*.a2g.reditools_table.txt.gz ); do 
            temp=${file%%.a2g*}; sample=${temp##*ls2/};
            tabix -T completeOutput/a2g_sites/every_a2g_site.bed $file | sed 's/$/\t'$sample'/g' >> completeOutput/a2g_sites/all_a2g.reditools_table.txt; 
        done

   awk '$0 ~ /AG/{print}' completeOutput/a2g_sites/all_a2g.reditools_table.txt > completeOutput/a2g_sites/allAG_a2g.reditools_table.txt
   
 
########## Split BAMs into F1R2 and F2R1 bams

mkdir completeOutput/merged_alignmentsSPLIT
module load samtools/1.9
for bamfile in completeOutput/merged_alignments/*.bam; do
	filename=$(basename $bamfile)
	experiment=${filename%%.a2g*}
	echo $filename; 
	mkdir completeOutput/merged_alignmentsSPLIT/${experiment}_negative
	mkdir completeOutput/merged_alignmentsSPLIT/${experiment}_positive
	samtools view -@ 8 -f 99 -b -o completeOutput/merged_alignmentsSPLIT/${experiment}_negative/${experiment}_negative99.bam ${bamfile};
	samtools view -@ 8 -f 147 -b -o completeOutput/merged_alignmentsSPLIT/${experiment}_negative/${experiment}_negative147.bam ${bamfile};
	
	samtools view -@ 8 -f 163 -b -o completeOutput/merged_alignmentsSPLIT/${experiment}_positive/${experiment}_positive163.bam ${bamfile};
	samtools view -@ 8 -f 83 -b -o completeOutput/merged_alignmentsSPLIT/${experiment}_positive/${experiment}_positive83.bam ${bamfile};
	
	samtools merge completeOutput/merged_alignmentsSPLIT/${experiment}_negative/${experiment}_negative.bam completeOutput/merged_alignmentsSPLIT/${experiment}_negative/*
	samtools merge completeOutput/merged_alignmentsSPLIT/${experiment}_positive/${experiment}_positive.bam completeOutput/merged_alignmentsSPLIT/${experiment}_positive/*
	
	#samtools depth -b all_edit_sites.bed ${bamfile} \
	#| sed 's/$/\t'${experiment}'/g' \
	#>> all_edit_sites_coveragedepth.tsv
done   
   

######find_hyperedited_regions

#xfindREDItoolPeaks script

#!/usr/bin/env bash
set -euo pipefail

# This script is designed to take in an output table from REDItools2 and detect "hyperedited regions", where 3 A-to-G changes fall within a $distance subsequence, and also collect any additional A-to-G changes that are less than $distance away from at least two other A-to-G changes in that peak.

substitution="${1:-AG}"
inputfile="${2:-/dev/stdin}"

<${inputfile} \
gawk -v FS="\t" -v OFS="\t" -v distance=20 -v mincount=5 -v substitution=${substitution} '
	# Skip one header line at the beginning.
	BEGIN { getline };
	# Select only lines with coverage above 5 reads and fraction A to G above 0.05.
	$5 > 4 && $8 == substitution && $9 > 0.05 {
		# If a transposition is close enough to the previous site, it will be added to the information about the current peak.
		if ( region == $1 && ($2 - secondlatest) < distance ) {
			secondlatest = latest;
			latest = $2;
			count++;
		# If a transposition is not close to the previous site, the current peak will be summarized and a new one started.
		# This will also trigger at the beginning of the file.
		} else {
			end = latest;
			if ( count >= mincount ) { print region, start, end, count };
			region = $1;
			start = $2;
			latest = $2;
			secondlatest = latest;
			count = 1;
		}
	}
	# Flush the buffer on every line for piping.
	{ fflush() }
	# Need to potentially output the final peak if it includes the very last AG.
	END {
		if ( count >= mincount ) { print region, start, end, count };
	}
'
#end

array=($(cat completeOutput/samples.txt))
for sample in $(seq 0 59); do zcat completeOutput/reditools2/${array[$sample]}.a2g.reditools_table.txt.gz | xfindREDItoolPeaks AG | gawk -v FS="\\t" -v OFS="\\t" '{{print $1, $2 - 1, $3, $1 ":" $2 "-" $3, $4 }} ' > completeOutput/hyperedited_peaks/${array[$sample]}.a2g.peaks.bed ; done
 
###### Feature counts

featureCounts -T 12 -t exon -g gene_name -a references/STAR/GRCh38/gencode.v41.primary_assembly.annotation.gtf -o genecounts.txt ./merged_alignments/*bam

 
 
