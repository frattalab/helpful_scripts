# helpful_scripts
 Helpful scripts for performing routine analyses


### check_rg_header.sh ###

Bash script that looks through directory and subdirectories for BAM files that are missing @RG flag at start of read group header and adds @RG accordingly. Used on some older processed data thathad these missing so re-processed BAMs can be used with more recent versions of samtools, pysam etc.
