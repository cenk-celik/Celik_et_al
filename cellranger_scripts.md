# Cell Ranger script for obtaining counts from the fastq files

- Cell Ranger works only on Linux. Make sure you have successfully installed the [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/6.1/installation). on your system.
- Place the *fastqs* in the same folder for each sample
- Download the reference genome for mouse from [10X Genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/6.1).

In shell, type the following to obtain counts for a sample:

```bash
cellranger count --id=sampleName \ # optional
--transcriptome=/home/users/ntu/cenkceli/yard/apps/refdata-gex-mm10-2020-A \ # reference genome directory
--fastqs=/scratch/users/ntu/cenkceli/sampleName_fastqs \ # the directory containing fastq files for a sample
--sample=sampleName-FLOWCELLID \ # must be the exact name of the fastqs before ***_S1***!
--expect-cells=8000 \ # how many cells were targeted
--include-introns=true # as Cell Ranger v6.1 and lesser do not include introns by default, add --include-introns=true if you need
```

If you have different FLOWCELLIDs for the same sample, grep all and feed into cellranger:

```bash
sample=sampleName
lsam=$(ls /scratch/users/ntu/cenkceli/sampleName_fastqs/ | grep $sample | awk -F['_'] '{print $1 "_" $2 "_" $3}'| uniq |xargs | sed -e 's/ /,/g')

cellranger count --id=$sample \
--fastqs=/scratch/users/ntu/cenkceli/sampleName_fastqs \
--sample=$lsam \
--transcriptome=/home/users/ntu/cenkceli/yard/apps/refdata-gex-mm10-2020-A \
--expect-cells=8000 \
--include-introns=true
```
