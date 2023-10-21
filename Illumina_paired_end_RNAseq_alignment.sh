#!/bin/bash

#module load python/3.9.8
#module load fastqc/0.11.9
#module load samtools/1.9
#module load trimgalore/0.6.7
#module load star/2.7.3a

workDir="xxxxx"

cd ${workDir}

mkdir ${workDir}/fastq_trim
mkdir ${workDir}/bam_xargs
mkdir ${workDir}/bamutil_xargs

ls ${workDir}/fastq | egrep R1_001.fastq.gz$ | xargs -P12 -I '{}' bash -c 'trim_galore --paired $2/fastq/$1 $2/fastq/${1%%_R1_001.fastq.gz}_R2_001.fastq.gz -o $2/fastq_trim' - '{}' ${workDir}


ls ${workDir}/fastq_trim | egrep R1_001_val_1.fq.gz$ | xargs -P2 -I '{}' bash -c 'STAR --runThreadN 8 \
	--genomeDir ~/Mouse/Gencode23_Star/Genome/ \
	--readFilesIn $2/fastq_trim/$1 $2/fastq_trim/${1%%_R1_001_val_1.fq.gz}_R2_001_val_2.fq.gz \
	--readFilesCommand zcat \
	--outSAMattributes All \
	--quantMode GeneCounts \
	--outSAMtype BAM SortedByCoordinate \
	--sjdbGTFfile ~/Mouse/Gencode23_Star/Genome/gencode.vM23.annotation.gtf \
	--chimSegmentMin 15 \
	--outSJfilterOverhangMin 15 15 15 15 \
	--alignSJoverhangMin 15 \
	--alignSJDBoverhangMin 15 \
	--seedSearchStartLmax 30 \
	--outFilterMultimapNmax 10 \
	--outFilterScoreMin 1 \
	--outFilterMatchNmin 1 \
	--outFilterMismatchNmax 2 \
	--chimScoreMin 15 \
	--chimScoreSeparation 10 \
	--chimJunctionOverhangMin 15 \
	--outWigType bedGraph \
	--outWigNorm None \
	--outFileNamePrefix $2/bam_xargs/${1%%_R[1|2]_001_val_[1|2].fq.gz}_' - '{}' ${workDir}


## index 

cd ${workDir}/bam_xargs

ls | egrep bam$ | xargs -P12 -I '{}' bash -c 'samtools index $1' - '{}'


### featureCounts 

cd ${workDir}
countsFileNam="counts_xargs.txt"

~/subread-2.0.0-Linux-x86_64/bin/featureCounts -T 12 \
-t exon \
-g gene_name \
-a ~/Mouse/Gencode23_Star/Genome/gencode.vM23.annotation.gtf \
-o ${workDir}/${countsFileNam} \
${workDir}/bam_xargs/*_Aligned.sortedByCoord.out.bam


## add ENSEMBL IDs

cd ${workDir}
countsFileNam="counts_xargs.txt"
bam_dir_reg="\/xxx\/xxx\/bam_xargs\/"


awk -v bam_dir=${bam_dir_reg} 'BEGIN{FS="\t"; OFS="\t"} FNR==NR && NR != 1 && length($2) > 0 {my_symbol[$2]=$2; my_ENSEMBL[$2]=$1; next} \
	FNR != NR && FNR == 1{next} \
	FNR != NR && FNR == 2{printf "%s\t%s\t", "symbol", "ENSEMBL"; \
		for (i=1; i<=6; i++){printf "%s\t", $i}; for (j=7; j< NF; j++){n=gensub(bam_dir, "", "g", $j); \
		printf "%s\t", n}; {printf "%s\n",  gensub(bam_dir, "", "g", $NF)}; next} \
		($1 in my_symbol){printf "%s\t%s\t", my_symbol[$1], my_ENSEMBL[$1]; \
		for (i=1; i < NF; i++){printf "%s\t", $i} printf "%s\n", $NF}' \
	~/mart/Mouse_mart_export/mart_export.txt ${countsFileNam} > counts_xargs_annot_ensembl.txt


### bamutilsj

cd ${workDir}/bam_xargs

ls | egrep bam$ | xargs -P12 -I '{}' bash -c '~/ngsutilsj bam-stats \
--gtf ~/Mouse/Gencode23_Star/Genome/gencode.vM23.annotation.gtf $1 2>&1 | tee $2/$1_bamutil.txt' - '{}' ${workDir}/bam_xargs




