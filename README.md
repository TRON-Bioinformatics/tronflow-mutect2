# TronFlow Mutect2

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/release/tron-bioinformatics/tronflow-mutect2?sort=semver)
[![Run tests](https://github.com/TRON-Bioinformatics/tronflow-mutect2/actions/workflows/automated_tests.yml/badge.svg?branch=master)](https://github.com/TRON-Bioinformatics/tronflow-mutect2/actions/workflows/automated_tests.yml)
[![DOI](https://zenodo.org/badge/355860788.svg)](https://zenodo.org/badge/latestdoi/355860788)
[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)
[![Powered by Nextflow](https://img.shields.io/badge/powered%20by-Nextflow-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://www.nextflow.io/)

The TronFlow Mutect2 pipeline is part of a collection of computational workflows for tumor-normal pair somatic variant calling.

Find the documentation here [![Documentation Status](https://readthedocs.org/projects/tronflow-docs/badge/?version=latest)](https://tronflow-docs.readthedocs.io/en/latest/?badge=latest)


This workflow implements the Mutect2 (Benjamin, 2019) best practices somatic variant calling of tumor-normal pairs.
![Mutect2 best practices](https://drive.google.com/uc?id=1rDDE0v_F2YCeXfQnS00w0MY3cAGQvfho)

It has the following steps:
* **Mutect2** - the somatic variant caller.
* **Learn read orientation model** - learn the prior probability of read orientation artifacts.
* **Pile-up summaries** - summarizes counts of reads that support reference, alternate and other alleles for given sites (optional).
* **Calculate contamination** - Given pileup data from GetPileupSummaries, calculates the fraction of reads coming from cross-sample contamination (optional).
* **Filter calls** - filters mutations from the raw Mutect2 variant calls
* **Funcotator annotation** - add functional annotations (optional)


## How to run it

Run it from GitHub as follows:
```
nextflow run tron-bioinformatics/tronflow-mutect2 -r v1.4.0 -profile conda --input_files $input --reference $reference --gnomad $gnomad
```

Otherwise download the project and run as follows:
```
nextflow main.nf -profile conda --input_files $input --reference $reference --gnomad $gnomad
```

Find the help as follows:
```
$ nextflow run tron-bioinformatics/tronflow-mutect2 --help

Usage:
    nextflow run tron-bioinformatics/tronflow-mutect2 -profile conda --input_files input_files [--reference reference.fasta]

This workflow is based on the implementation at /code/iCaM/scripts/mutect2_ID.sh

Input:
    * input_files: the path to a tab-separated values file containing in each row the sample name, tumor bam and normal bam
    The input file does not have header!
    Example input file:
    name1	tumor_bam1	normal_bam1
    name2	tumor_bam2	normal_bam2
    * reference: path to the FASTA genome reference (indexes expected *.fai, *.dict)
    
Optional input:
    * input_name: sample name (alternative to --input_files)
    * input_tumor_bam: comma separated list of tumor BAMs (alternative to --input_files)
    * input_normal_bam: comma separated list of normal BAMs (alternative to --input_files)
    * gnomad: path to the gnomad VCF or other germline resource (recommended). If not provided the contamination will 
    not be estimated and the filter of common germline variants will be disabled
    * pon: path to the panel of normals VCF
    * intervals: path to a BED file containing the regions to analyse
    * output: the folder where to publish output (default: output)
    * enable_bam_output: outputs a new BAM file with the Mutect2 reassembly of reads (default: false)
    * disable_common_germline_filter: disable the use of GnomAD to filter out common variants in the population
    from the somatic calls. The GnomAD can still be provided though as this common SNPs are used elsewhere to
    calculate the contamination (default: false)
    * funcotator: To use Funcotator, supply the path to a database to be used. (can be downloaded from GATK FTP server)
    * reference_version_funcotator: version of the reference genome (default: "hg19")
    * output_format_funcotator: the output format of Funcotator. Can be VCF or MAF (default: "MAF")
    * transcript_selection_mode_funcotator: transcript selection method can be CANONICAL, BEST_EFFECT or ALL. (default: CANONICAL)
    * memory_mutect2: the ammount of memory used by mutect2 (default: 16g)
    * memory_read_orientation: the ammount of memory used by learn read orientation (default: 16g)
    * memory_pileup: the ammount of memory used by pileup (default: 32g)
    * memory_contamination: the ammount of memory used by contamination (default: 16g)
    * memory_filter: the ammount of memory used by filter (default: 16g)
    * memory_funcotator: the ammount of memory used by filter (default: 16g)
    * args_filter: optional arguments to the FilterMutectCalls function of GATK (e.g.: "--contamination-estimate 0.4 --min-allele-fraction 0.05 --min-reads-per-strand 1 --unique-alt-read-count 4") (see FilterMutectCalls documentation)
    * args_funcotator: optional arguments to Funcotator (e.g. "--remove-filtered-variants true")  (see Funcotator documentation)
    * args_mutect2: optional arguments to Mutect2 (e.g. "--sites-only-vcf-output")  (see Mutect2 documentation)

Output:
    * Output VCF
    * Other intermediate files
```


### Input tables

The table with BAM files expects three tab-separated columns without a header.
Multiple tumor or normal BAMs can be provided separated by commas.

| Sample name          | Tumor BAMs                      | Normal BAMs                  |
|----------------------|---------------------------------|------------------------------|
| sample_1             | /path/to/sample_1_tumor.bam      |    /path/to/sample_1_normal.bam   |
| sample_2             | /path/to/sample_2_tumor_1.bam,/path/to/sample_2_tumor_2.bam      |    /path/to/sample_2_normal.bam,/path/to/sample_2_normal_2.bam   |

### About read group tags in BAM headers

Mutect2 relies on several read group tags to be present in the BAM header. 
The commpulsory tags are: `RG:ID`, `RG:PU`, `RG:SM`, `RG:PL` and `RG:LB`.
If your BAM files do not have these read group tags use 
[Picard's AddOrReplaceReadGroups](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-).

There are some further constraints in the sample tag (`RG:SM`) to distinguish normal and tumor samples.
Hence, this workflow expects that tumor and normal BAMs have different values of RGSM; 
and when replicates are provided all normal BAMs must have the same RGSM; and the same applies for all tumor BAMs.
The workflow will fail if these constraints are not met.


## Resources

- FASTA reference genome with fai and dict indexes (see https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format for instructions on building the indices)
- Analysis intervals file in BED format. These intervals determine the regions where variants will be called
- VCF file with common germline variants (see [GnomAD](#gnomad))
- Optionally, a panel of normals (PON) may be used (see [PON](#pon))

### GnomAD

GnomAD (Karczewski, 2020) is the standard de facto database for germline variants population allele frequencies. Mutect2 employs GnomAD as prior knowledge to reject potential germline variants and it also uses the SNP to estimate contamination.

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

### Panel Of Normals (PON)

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

Once the panel of normals is created pass it to the workflow using the parameter `--pon`.

### Configuring Funcotator

Funcotator annotation is an optional step. To configure funcotator follow the indications here https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial.

In order to use funcotator provide the path to your local funcotator database with the parameter `--funcotator`.
Also, make sure that the reference version provided to funcotator with `--reference_version_funcotator` is consistent with the provided reference with `--reference`. 


## How to run the Panel of Normals (PON) pipeline

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

## References

- Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316–319. https://doi.org/10.1038/nbt.3820
- Benjamin, D., Sato, T., Cibulskis, K., Getz, G., Stewart, C., & Lichtenstein, L. (2019). Calling Somatic SNVs and Indels with Mutect2. BioRxiv. https://doi.org/10.1101/861054
- GATK team. Somatic short variant discovery (SNVs + Indels). Retrieved from https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-
- Karczewski, K.J., Francioli, L.C., Tiao, G. et al. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature 581, 434–443 (2020). https://doi.org/10.1038/s41586-020-2308-7
