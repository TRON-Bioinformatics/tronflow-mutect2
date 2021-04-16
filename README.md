# TronFlow Mutect2

A nextflow pipeline implementing Mutect2 best practices somatic variant calling of tumor-normal pairs. See https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-

![Mutect2 best practices](mutect2_best_practices.png)


## How to run it

```
$ nextflow run tron-bioinformatics/tronflow-mutect2 -profile conda --help
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
    * memory: the ammount of memory used by each job (default: 16g)
    * cpus: the number of CPUs used by each job (default: 2)
    * disable_common_germline_filter: disable the use of GnomAD to filter out common variants in the population
    from the somatic calls. The GnomAD resource is still required though as this common SNPs are used elsewhere to
    calculate the contamination (default: false)

Output:
    * Output VCF
    * Other intermediate files
```

## Requirements

- GATK 4.2.0.0
- Java 8

All dependencies are set in a conda environment, thus no installation is required if conda is used.

## Resources and data requirements

- FASTA reference genome with fai and dict indexes (see https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format for instructions on building the indices)
- Analysis intervals file in BED format. These intervals determine the regions where varianst will be called
- VCF file with common germline variants (see [GnomAD](#gnomad))
- Input BAM files require that read groups are added to them, furthermore tumor and normal must have different read groups (see [Adding read groups to BAM files](#adding-read-groups-to-bam-files)).
- Optionally, a panel of normals (PON) may be used (see [PON](#pon))

### GnomAD

GnomAD is the standard de facto database for germline variants population allele frequencies. Mutect2 employs GnomAD as prior knowledge to reject potential germline variants and it also uses the SNP to estimate contamination.

If you want to disable the use of GnomAD to filter out common germline variants in the somatic calls use `--disable_common_germline_filter`, GnomAD will still be used to estimate the contamination.

This resource has a total of 14,967,411 variants, of which 14,078,157 SNVs and 889,254 indels. No variants reported in mitochondrial chromosome. Frequencies are mostly low as expected, although there are some variants with a frequency of 1.0. Overall, we have 95,542 common SNVs (ie: AF > 5%), 96,599 low frequency SNVs (ie: AF<=5% and AF >= 0.5%) and 13,886,016 rare SNVs (AF < 0.5%) (of which 13,703,545 have AF < 0.1%); and 13,200 common indels, 12,444 low frequency indels and 863,610 rare indels.

GnomAD v2.1 for the coding region can be downloaded from https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz

The details on this release are described here https://macarthurlab.org/2018/10/17/gnomad-v2-1/

We keep only the variants passing all filters. We remove all annotations except AC, AF and AN, as Mutect does not use any other annotations such as specific population frequencies.
```
bcftools annotate --include 'FILTER="PASS"' --remove ^INFO/AC,INFO/AF,INFO/AN /projects/data/gatk_bundle/b37/gnomad.exomes.r2.1.1.sites.vcf.bgz --output-type z --output /projects/data/gatk_bundle/b37/gnomad.exomes.r2.1.1.sites.PASS.only_af.vcf.bgz --threads 4
```

GnomAD file is in b37, thus it may be needed to lift over to hg19 for instance. Lift over chain files can be found here: ftp://gsapubftp-anonymous@ftp.broadinstitute.org/Liftover_Chain_Files.
```
java -jar /code/picard/2.21.2/picard.jar LiftoverVcf INPUT=/projects/data/gatk_bundle/b37/gnomad.exomes.r2.1.1.sites.PASS.only_af.vcf.bgz OUTPUT=/projects/data/gatk_bundle/hg19/gnomad.exomes.
```

### Adding read groups to BAM files

Use Picard's AddOrReplaceReadGroups tool (see https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-).
You will need to set a different sample name for tumor and normal (parameter `-SM`), you can just use `tumor` and `normal` in absence of a better naming.

### PON

The PON is used to filter out technical artifacts from the somatic variant calls. The PON is recommended to be formed from technically similar samples (ie: same sequencing platform, same sample preparation), from healthy and young individuals and to be formed by a minimum of 40 samples (see https://gatkforums.broadinstitute.org/gatk/discussion/11053/panel-of-normals-pon ).

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

This is implemented in the pipeline `mutect2_pon.vcf`.


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
