Helpful scripts for performing routine analyses


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

Fix collect_input_stat to operate more flexibly
Mapping_stat.py should read from the command line not global variables
