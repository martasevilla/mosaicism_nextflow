# ![CIBERER/GdTBioinfo-nf-mosaicism](docs/images/ciberer_logo.png)

<!--[![GitHub Actions CI Status](https://github.com/CIBERER/GdTBioinfo-nf-mosaicism /workflows/nf-core%20CI/badge.svg)](https://github.com/CIBERER/GdTBioinfo-nf-mosaicism /actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/CIBERER/GdTBioinfo-nf-mosaicism /workflows/nf-core%20linting/badge.svg)](https://github.com/CIBERER/GdTBioinfo-nf-mosaicism /actions?query=workflow%3A%22nf-core+linting%22) -->

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/CIBERER/GdTBioinfo-nf-mosaicism )

<!-- [![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23mosaicism-4A154B?logo=slack)](https://nfcore.slack.com/channels/mosaicism) -->
[![Follow on Twitter](https://img.shields.io/badge/twitter-%40CIBERER-1DA1F2?logo=twitter)](https://twitter.com/CIBERER)
[![Watch on YouTube](https://img.shields.io/badge/youtube-CentrodeInvestigaci%C3%B3nCIBER-FF0000?logo=youtube)](https://www.youtube.com/c/CentrodeInvestigaci%C3%B3nCIBER)

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

**CIBERER/GdTBioinfo-nf-mosaicism** is a bioinformatics best-practice analysis pipeline for detection of mosaicism through NGS data.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

This pipeline uses two variant callers to detect low frequency variants, VarScan and VarDict. The result file is a VCF file obtained after intersecting the two VCF files.

1. VarScan ([`VarScan`](https://varscan.sourceforge.net/))

1.1. Sorting BAM file ([`Samtools`](http://www.htslib.org/))

1.2. Obtaining mpileup file ([`Samtools`](http://www.htslib.org/))

1.3. Variant Calling ([`VarScan`](https://varscan.sourceforge.net/))

2. VarDict ([`VarDictJava`](https://github.com/AstraZeneca-NGS/VarDictJava))

2.1. Variant Calling ([`VarDictJava`](https://github.com/AstraZeneca-NGS/VarDictJava))

3. Mutect2 ([Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/21905083931035-Mutect2))

3.1 Generate reference sequence dictionary ([CreateSequenceDictionary](https://gatk.broadinstitute.org/hc/en-us/articles/21905059452187-CreateSequenceDictionary-Picard))

3.2 Variant Calling ([Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/21905083931035-Mutect2))

4. Intersection of the three VCF files ([`Bedtools`](https://bedtools.readthedocs.io/en/latest/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow pull CIBERER/GdTBioinfo-nf-mosaicism
   nextflow run CIBERER/GdTBioinfo-nf-mosaicism -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Then, download the [profile_test.md5](./tests/profile_test.md5) file, locate it in the root folder and run the next command:

   ```console
   md5sum -c profile_test.md5
   ```

   If the output is like this, everything is working fine!:

   ```console
   ./results/join/Sample_zenodo_T1.vcf: OK
   ./results/mutect2/Sample_zenodo_T1.mutect2.vcf.gz: OK
   ./results/mutect2/Sample_zenodo_T1.mutect2.vcf.gz.stats: OK
   ./results/mutect2/Sample_zenodo_T1.mutect2.vcf.gz.tbi: OK
   ./results/vardictjava/Sample_zenodo_T1.mpileup: OK
   ./results/vardictjava/Sample_zenodo_T1.txt.gz: OK
   ./results/vardictjava/Sample_zenodo_T1..vcf.gz: OK
   ./results/varscan/Sample_zenodo_T1.v2.vcf.gz: OK
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```console
   nextflow run CIBERER/GdTBioinfo-nf-mosaicism --input path/to/samplesheet.csv --outdir <OUTDIR> --fasta path/to/reference/genome.fa --fai path/to/reference/genome.fa.fai --chrom_sizes path/to/referece/chrom.sizes --germline_resource path/to/population/germline_resource.vcf.gz --germline_resource_tbi path/to/population/germline_resource.vcf.gz.tbi -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Documentation

The  CIBERER-pipelines /mosaicism pipeline comes with documentation about the pipeline [usage](./docs/usage.md), [parameters](./schema.md) and [output](./docs/output.md).


## Credits

This pipeline was originally written by Carlos Ruiz, Marta Sevilla and Yolanda Benítez at the [Centro de Investigaciones Biomédicas en Red de Enfermedades Raras, CIBERER](https://www.ciberer.es/).

Main authors:

- [Carlos Ruiz](https://github.com/yocra3)
- [Marta Sevilla](https://github.com/martasevilla)
- [Yolanda Benítez](https://github.com/yolandabq)
- [Pedro Garrido](https://github.com/pedro-garridor)
- [Paula Iborra](https://github.com/paulaidt)

## Contributions and Support

If you would like to contribute to this pipeline, please contact with:

carlos.ruiza@upf.edu, marta.sevilla@upf.edu or yolanda.benitez@ciberer.es

 <!-- For further information or help, don't hesitate to get in touch on the [Slack `#mosaicism` channel](https://nfcore.slack.com/channels/mosaicism) (you can join with [this invite](https://nf-co.re/join/slack)).
-->

## Citations

1. VarScan ([`VarScan`](https://varscan.sourceforge.net/))
2. Samtools ([`Samtools`](http://www.htslib.org/))
3. VarDictJava ([`VarDictJava`](https://github.com/AstraZeneca-NGS/VarDictJava))
4. Mutect2 ([GATK](https://gatk.broadinstitute.org/hc/en-us))
5. Bedtools [`Bedtools`](https://bedtools.readthedocs.io/en/latest/)
6. Tabix([`Tabix`](http://www.htslib.org/doc/tabix.html))


<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  CIBERER/GdTBioinfo-nf-mosaicism  for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
 
> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
> In addition, references of tools and data used in this pipeline are as follows:
