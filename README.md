# chromo-vis

R script for plotting variants, coverage depth, and runs of homozygozity on individual chromosomes.

![Example image](https://raw.githubusercontent.com/MNMdiagnostics/chromo-vis/master/example.png)

### Prerequisites:
 
The script requires following R packages:
  - VariantAnnotation
  - BSgenome (with hg38 installed)
  - optparse

### Usage:
  
  To plot variants, a tabixed VCF file with AD and DP fields is required. 
  Only variants with PASS value in the FILTER column are used.
  Currently, only hg38 coordinates system is supported, but create an issue if you need hg19 too.

  Running
  
  `./chromovis.R -v mysample.vcf.gz`
  
  generates 24 images (one for each of chr1-22, X and Y) named `mysample.vcf.gz.chr[1-22XY].png`. 
  
  Images can be saved with a custom prefix (`-o`), contain regions of homozygozity from supplied BED file (`-r`), and be generated for selected chromosomes (`-c`), e.g.
  
  `./chromovis.R -v mysample.vcf.gz -r mysample.roh.bed -o chromovis/mysample -c chr5,chr6,chrX`
  
  To list all available options run
  
  `./chromovis.R --h`
  
  


