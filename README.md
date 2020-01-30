# helpful_scripts
 Helpful scripts for performing routine analyses


### check_rg_header.sh ###

Bash script that looks through directory and subdirectories for BAM files that are missing @RG flag at start of read group header and adds @RG accordingly. Used on some older processed data thathad these missing so re-processed BAMs can be used with more recent versions of samtools, pysam etc.

### bed2gtf.py ###

Takes a basic BED file and converts it into a GTF format which can be used by featureCounts

### collect_input_stat.sh ###

Takes a folder of aligned bams (right now pretty inflexibly) and runs the picard tools collect insertsize metrics tools on them

### collect_input_stat.sh ###
