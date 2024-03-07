#!/bin/bash

samtools fastq \
  -0 /dev/null \
  -1 NA12878_1.fastq.gz \
  -2 NA12878_2.fastq.gz \
  -s /dev/null \
  mosaicism_nextflow/testdata/my.bam
bwa mem \
  -t 12 \
  -M \
  -R "@RG\tID:NA12878\tSM:NA12878\tPL:illumina\tLB:NA12878" \
  chr17_var_regions.fa \
  NA12878_1.fastq.gz NA12878_2.fastq.gz | \
samtools sort --no-PG -O bam -o sorted.bam
samtools index sorted.bam
java -jar ~/picard.jar MarkDuplicates \
  I=sorted.bam \
  O=my.bam \
  M=/dev/null
samtools index my.bam
rm sorted.bam*
mv my.bam* mosaicism_nextflow/testdata/