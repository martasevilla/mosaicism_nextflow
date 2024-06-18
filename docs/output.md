# CIBERER-pipelines/mosaicism: Output

# Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

# Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Samtools](#samtools) - It prepares the bam file for VarDictJava. Its final output is a zipped mpileup file.   
- [Tabix](#tabix) - It unzips the mpileup file.
- [VarDictJava](#vardictjava) - VCF file obtained by VarDict variant caller.
- [VarScan](#varscan) - VCF file obtained by VarScan variant caller.
- [CreateSequenceDictionary](#createsequencedictionary) - Generate genome dictionary for Mutect2.
- [Mutect2](#mutect2) - VCF file obtained by Mutect2 variant caller.
- [Bedtools](#bedtools) - Final VCF file obtained after intersecting both VCF files.  
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

# Steps and output files
<!---
### Samtools

<details markdown="1">
<summary>Output files</summary>

- `Samtools/`
  - `*.bam`: Sorted bam file.
  - `*.mpileip.gz`: Zip archive containing the mpileup file (VarDictJava input).

</details>

[Samtools](https://varscan.sourceforge.net/) prepares the data input for VarDictJava. This variant calling needs a mpileup file. It produces a zipped mpileip file, that is unzipped by [Tabix](http://www.htslib.org/doc/tabix.html). For further reading and documentation see the [Samtools help pages](http://www.htslib.org/doc/#manual-pages).
-->

### VarDictJava

<details markdown="1">
<summary>Input files</summary>

- `bam file`
- `bam.bai file`
- `bed file containing the regions of interest`
- `reference fasta file`
- `reference fasta.fai file`

</details>

<details markdown="1">
<summary>Output files</summary>

- `vardictjava/*.vcf.gz`

</details>

[VarDictJava](https://github.com/AstraZeneca-NGS/VarDictJava) is an ultra sensitive variant caller for both single and paired sample variant calling from BAM files.

### VarScan

<details markdown="1">
<summary>Input files</summary>

- `mpileup file`

</details>

<details markdown="1">
<summary>Output files</summary>

- `varscan/*.vcf.gz`

</details>

[VarScan](https://varscan.sourceforge.net/) is a platform-independent mutation caller that can be used to detect germline variants. It employs a robust heuristic/statistic approach to call variants that meet desired thresholds for read depth, base quality, variant allele frequency, and statistical significance.

The parameters used in VarScan can be found in the module file of [VarScan](../modules/local/varscan/main.nf).

```console
varscan mpileup2snp ${mpileup} --min-var-freq 0.05 --p-value 1 --output-vcf 1
```

<!--
### CreateSequenceDictionary

<details markdown="1">
<summary>Input files</summary>

- `reference fasta file`

</details>

<details markdown="1">
<summary>Output files</summary>

- `CreateSequenceDictionary/*.dict`

</details>

[CreateSequenceDictionary](https://gatk.broadinstitute.org/hc/en-us/articles/21905059452187-CreateSequenceDictionary-Picard) creates a sequence dictionary for a reference sequence. This tool creates a sequence dictionary file (with ".dict" extension) from a reference sequence provided in FASTA format, which is required by many processing and analysis tools. The output file contains a header but no SAMRecords, and the header contains only sequence records.
-->

### Mutect2

<details markdown="1">
<summary>Input files</summary>

- `bam file`
- `bam.bai file`
- `bed file containing the regions of interest`
- `reference fasta file`
- `reference fasta.fai file`
- `reference dict file`
- `vcf germline resouce file`
- `vcf.tbi germline resource file`

</details>

<details markdown="1">
<summary>Output files</summary>

- `mutect2/`
  - `*.vcf.gz`
  - `*.vcf.gz.tbi`
  - `*.vcf.gz.stats`

</details>

[Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/21905083931035-Mutect2) calls somatic short mutations via local assembly of haplotypes. Short mutations include single nucleotide (SNA) and insertion and deletion (indel) alterations. The caller uses a Bayesian somatic genotyping model that differs from the original MuTect by [Cibulskis et al., 2013](http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html) and uses the assembly-based machinery of [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php).

### Bedtools

<details markdown="1">
<summary>Input files</summary>

- `vcf file 1`
- `vcf file 2`
- `vcf file 3`

</details>

<details markdown="1">
<summary>Output files</summary>

- `join/*.vcf.gz` : vcf containing the intersection between the variant callers.

</details>

[Bedtools](https://bedtools.readthedocs.io/en/latest/) perform the intersection of the vcf obtained by the three variant callers. This way, we keep the variants that are more likely to be real variants.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
