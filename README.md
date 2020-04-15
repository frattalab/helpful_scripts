Helpful scripts for performing routine analyses

### feature-annotation-overlap.R ###

Basically bedtools intersect but in R. Script that can take an annotation file (GTF/BED, or can use inbuilt annotatr genic annotations) and BED file containing coordinates of features (e.g. iCLIP sites) and report where feature coordinates intersect with annotation coordinates. Outputs tab-separated tablewhere each line contains feature info and info of annotation which which it overlaps. Can also filter for overlaps by 'feature field' (3rd column/field) in provided GTF file (e.g. gene,exon,intron,UTR etc.).

Run ```Rscript feature-annotation-overlap.R -h``` at the command line to get more info on options.

(Stay tuned might add more functionality)

### check_rg_header.sh ###

Bash script that looks through directory and subdirectories for BAM files that are missing @RG flag at start of read group header and adds @RG accordingly. Used on some older processed data thathad these missing so re-processed BAMs can be used with more recent versions of samtools, pysam etc.

### bed2gtf.py ###

Takes a basic BED file and converts it into a GTF format which can be used by featureCounts

### collect_input_stat.sh ###

Takes a folder of aligned bams (right now pretty inflexibly) and runs the picard tools collect insertsize metrics tools on them

### make_deseq_dfs.R ###

Take a data table or data frame and returns a list of a data frame with rownames as genes, and a metadata table for use with the DESEQ2 R tool

### create_feature_countable.R ###

Takes a folder path with featureCount tables per sample and makes it into a tidy table

### converting_gene_lists.R ###

Helper R functions for converting between mouse/human(Using bioMart) and converting between ENTREZ,HUGO,ENSEMBL ids using org.Annotation.dbs

### librarySize.py ###

Functions which read the STARLogFiles and FeatureCountsMapped files to retrieve mapping statistics

### mapping_stats.py ###

Imports the functions from librarySize and creates a csv with mapping stats based on the variables in the script. Currently a kludge solution

### TODO's ###

Streamline feature-annotation-overlap.R command line options
Test feature-annotation-overlap.R with -d option to use inbuilt annotations
Fix collect_input_stat to operate more flexibly
Mapping_stat.py should read from the command line not global variables
