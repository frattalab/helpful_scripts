#annotatr stuff from Anna-Leigh Brown, extra 'functionality' by Sam Bryce-Smith

knitr::opts_chunk$set(echo = TRUE)

if(!require("optparse")) install.packages("optparse")
library(optparse)

option_list = list(
  make_option(c("-d","--download"), type = "character",
              default = NULL,
              help = "Include if want to automatically download GTF/gene annoations instead of providing. Must be followed by one of c(hg19, hg38, mm9, mm10, rn4, rn5, rn6, dm3, and dm6) - these are built-in genic annotations in annotatr package. Exclude option if providing own annotation)"),
  make_option(c("-a","--annotation"), type = "character",
              default = NULL,
              help = "follow with path to annotation file containing genomic features of interest (BED,GTF). Don't include option if want to script to download annotation file"
              ),
  make_option(c("--gtf"),
              action= "store_true",
              help = "Specify if provided annotation file is in GTF/GFF format. Do not include if using '-d/--download' options or input annotation is BED format"),
  make_option(c("--bed"),
              action = "store_true",
              help = "Specify if provided annotation file is BED format. Do not include if using '-d/--download' options or input annotation is GTF format (-gtf)"),
  make_option(c("-f","--features"),
              default = NULL,
              help = "follow with path to BED file containing iCLIP coordinates/feature coordinates with which to intersect the annotation file"),
  make_option(c("-s", "--strand-specific"),
              action = "store_true",
              help = "specify whether want intersect to be strand-specific (feautures and annotations will only be reported as overlapping if on the same strand. Don't include if want to perform non-strand-specific overlap."),
  make_option(c("--filter"), type = "character",
              default = NULL,
              help = "optional parameter to filter overlap based on feature type e.g. c(gene,exon,intron,UTR) etc. (GTF annotation only). Make sure matches with name present in 'feature' field of GTF/GFF file (3rd column/field)"),
  make_option(c("-o","--output"), type = "character",
              default = NULL,
              help = "Optional - specify path to output directory - MUST END WITH '/'. Do not include if want to output to current working directory. [default=%default]"))


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

##If option is -h (or --help) print help messages and exit
if(isTRUE(opt$help)){
  print_help(opt_parser)
  stop(call. = FALSE)
}

##Make sure provided options are compatible
#1. Download option + annotation, gtf & bed not set (can't remember why using this one...)
if(isFALSE(opt$download) && is.null(opt$annotation) && isFALSE(opt$gtf) && isFALSE(opt$bed)) {
  stop("-d option & c(annotation,gtf,bed) options are mutually exclusive. Only provide -d if want to use inbuilt annotations. Do not include -d option if want to provide annotations")
}

#2 annotation file provided and bed or GTF option not provided
if(!(is.null(opt$annotation)) && isFALSE(opt$gtf) && isFALSE(opt$bed)) {
  
}

#Print what are input files
print(paste("annotation file is ",opt$annotation,sep = ""))
print(paste("feature file is ",opt$features, sep = ""))






#Script to test for overlap with between genomic regions of interest and iCLIP-derived protein binding sites
#Can further test for enrichment of binding in a provided gene list

#Required input:

#GTF file containing genomic coordinates of interest with which want to overlap
#BED file containing coordinates of processed iCLIP clusters (of protein of interest)

if (!require("pacman")) install.packages("pacman")
library(pacman)
#library(GenomicRanges)

pacman::p_load(annotatr,GenomicRanges,rtracklayer)

#Load in annotation file (features e.g. gene coordinates etc.) to make a GRanges object
if (isTRUE(opt$gtf)) {
  annotation = import(opt$annotation)
} else if(opt$bed) {
  annotation = annotatr::read_regions(opt$annotation)
} else if(!is.null(opt$download)) {
  annotation = build_annotations(genome = opt$download, annotations = paste(opt$download,"_basicgenes", sep = ""))
}

#Read in provided BED file containing coordinates of features wish to intersect with annotation
feature_peaks <- read_regions(opt$features)

#Perform overlap/intersect - creates GRanges object reporting every intersect between iCLIP & gtf
if(isTRUE(opt$`strand-specific`)) {
  annotation_feature_overlap <- annotate_regions(regions = GRanges(feature_peaks), annotations = annotation, ignore.strand = FALSE)
} else if(isFALSE(opt$`strand-specific`)) {
  annotation_feature_overlap <- annotate_regions(regions = GRanges(feature_peaks), annotations = annotation, ignore.strand = TRUE)
}

#Convert Granges to dataframe - can manipulate more easily (and write to file)
df_annotation_feature_overlap <- as.data.frame(annotation_feature_overlap)

#Filter overlaps by feature (if set)
if(!is.null(opt$filter)) {
  if(!require("dplyr")) install.packages("dplyr")
  library(dplyr)
  #annotate_regions creates table with intersecting regions between feature file & annotation file - reports feature columns first followed by annotation columns
  #annotation columns named 'annot.<field_name>' e.g. 'annot.strand'
  #This case - 'annot.type' is feature column name we're interested in
  #"annot.type"
  
  #only do filtering if can find provided feature type string in annotation column
  if(any((opt$filter %in% df_annotation_feature_overlap[,"annot.type"]),TRUE)) {
    df_filtered_annotation_feature_overlap = df_annotation_feature_overlap %>%
      filter(annot.type == opt$filter)
  }
  else {
    print(paste("could not find",opt$filter," in the annot.type column. Check if provided string is present in the provided GTF file", sep = ""))
  }
  
}


#write intersect table to file (tab-delimited)
write.table(df_annotation_feature_overlap, file = paste(opt$output,"annotation_feature_overlap.tsv", sep = ""), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

#Write filtered table to file if created (tab-delimited)
if(exists("df_filtered_annotation_feature_overlap")) {
  write.table(df_filtered_annotation_feature_overlap, file = paste(opt$output,"filtered_annotation_feature_overlap.tsv",sep = ""), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
}

#If used custom annotation - print the built object to file
if(!is.null(opt$download)) {
  write.table(annotation, file = paste(opt$output,"annotation.tsv", sep = ""), row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
}

