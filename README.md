# sppIDer_TGS
*sppIDer* for *T*hird *G*eneration *S*equencing data
This is a updated version of original [sppIDer](https://github.com/GLBRC/sppIDer).


## what is new?

### 1. adapted mapping tool to third generation sequencing data
sppIDer using bwa to map the read to combined reference, with not suitable for third generation sequencing data.
So we add new option to it:
if you still want use bwa mem to mapping your PacBio / Nanopore data to reference, you can use

```
new_sppIDer.py \
  --out OUT --ref REF \
  --r1 oacbio.fastq \
  --seq-type PacBio \
  --mapping-tool bwa \
  --cores 12
```
the bwa will run with `bwa mem -x pacbio` parameter

### 2. add new function 
- new you can set how many core used in your analysis.
- if you want using sppIDer to classify your hybrid speceis data, you can use --seq-name option to output it.(still under working)
- more todo


## TODO
1. the output fig are not publishable. We will update the code to improve it.
2. we will provide new classify strategies, to improve the performance and new functions!
3. add download function, you can download reference data from NCBI easy and fast!
4. add sppIDer_TGS to bioconda channel for easy installation.