# tronflow-mutect2

A nextflow pipeline for Mutect2 best practices somatic variant calling of tumor-normal pairs.

## Requirements and resources

Latest version is integrated with GATK and the latest release is 4.1. It requires Java 8.

Data has similar requirements to GATK and to Mutect 1, such as having read groups set in the BAM. But, given that Mutect 2 is based on haplotype assembly the realignment around indels is not required in BAM preprocessing.

### GnomAD

GnomAD is the standard de facto database for germline variants population allele frequencies. Mutect2 employs GnomAD as prior knowledge to reject potential germline variants.

This resource has a total of 14,967,411 variants, of which 14,078,157 SNVs and 889,254 indels. No variants reported in mitochondrial chromosome. Frequencies are mostly low as expected, although there are some variants with a frequency of 1.0. Overall, we have 95,542 common SNVs (ie: AF > 5%), 96,599 low frequency SNVs (ie: AF<=5% and AF >= 0.5%) and 13,886,016 rare SNVs (AF < 0.5%) (of which 13,703,545 have AF < 0.1%); and 13,200 common indels, 12,444 low frequency indels and 863,610 rare indels.

GnomAD v2.1 for the coding region is downloaded from https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz

The details on this realease are described here https://macarthurlab.org/2018/10/17/gnomad-v2-1/

We keep only the variants passing all filters. We remove all annotations except AC, AF and AN, as Mutect does not use any other annotations such as specific population frequencies.
```
bcftools annotate --include 'FILTER="PASS"' --remove ^INFO/AC,INFO/AF,INFO/AN /projects/data/gatk_bundle/b37/gnomad.exomes.r2.1.1.sites.vcf.bgz --output-type z --output /projects/data/gatk_bundle/b37/gnomad.exomes.r2.1.1.sites.PASS.only_af.vcf.bgz --threads 4
```

GnomAD file in b37 will be lifted over to hg19, which implies just changing the chromosome names as GnomAD only provides data on the canonical and sexual chromosomes. Lift over chain file from b37 to hg19: ftp://gsapubftp-anonymous@ftp.broadinstitute.org/Liftover_Chain_Files/b37tohg19.chain. All variants were lifted over.
```
java -jar /code/picard/2.21.2/picard.jar LiftoverVcf INPUT=/projects/data/gatk_bundle/b37/gnomad.exomes.r2.1.1.sites.PASS.only_af.vcf.bgz OUTPUT=/projects/data/gatk_bundle/hg19/gnomad.exomes.
```

### Panel of normals (PON)

The PON is used to filter out germline-like variants from the somatic variant calls. The PON is recommended to be formed from technically similar samples (ie: same sequencing platform, same sample preparation), from healthy and young individuals and to be formed by a minimum of 40 samples (see https://gatkforums.broadinstitute.org/gatk/discussion/11053/panel-of-normals-pon ).

The normal samples are processed by the BAM preprocessing pipeline including marking duplicates and BQSR.

Run MuTect2 on each normal sample as follows:

```
java -Xmx16g -jar /code/gatk/4.1.3.0/gatk-package-4.1.3.0-local.jar \
    Mutect2 \
    --reference ${params.reference} \
    --intervals ${interval} \
    --input ${bam} \
    --tumor-sample ${name} \
    --max-mnp-distance 0 \
    --output ${bam.baseName}.${interval.baseName}.mutect.vcf
```

Note the parameter "--max-mnp-distance 0" is needed to avoid MNPs being called.

The multiple VCFs need to be combined with the GATK tool "CreateSomaticPanelOfNormals".

This is implemented in the pipeline mutect2_pon.vcf.


## Best practices

This pipeline implements the best practices described here: https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146 and https://software.broadinstitute.org/gatk/documentation/article?id=24057

![Mutect2 best practices](mutect2_best_practices.png)


## How to run it

```
$ nextflow main.nf --help
Usage:
    nextflow main.nf --input_files input_files

This workflow is based on the implementation at /code/iCaM/scripts/mutect2_ID.sh

Input:
    * input_files: the path to a tab-separated values file containing in each row the sample name, tumor bam and normal bam
    The input file does not have header!
    Example input file:
    name1	tumor_bam1	normal_bam1
    name2	tumor_bam2	normal_bam2

Optional input:
    * reference: path to the FASTA genome reference (indexes expected *.fai, *.dict)
    * intervals: path to a BED file containing the regions to analyse
    * gnomad: path to the gnomad VCF
    * NOTE: if any of the above parameters is not provided, default hg19 resources will be used
    * output: the folder where to publish output

Output:
    * Output VCF
```

## How to run the PON pipeline

```
$ nextflow mutect2_pon.nf --help
Usage:
    mutect2_pon.nf --input_files input_files

This workflow aims to compute a panel of normals to be used with MuTect2

Input:
    * input_files: the path to a file containing in each row the sample name and the path to a BAM file to be included in the PON
    	example:
	sample1	/path/to/sample1.bam
	sample2	/path/to/sample2.bam
	NOTE: the sample name must be set in the @SN annotation

Optional input:
    * output: the folder where to publish output

Output:
    * Output combined VCF pon.vcf
```
