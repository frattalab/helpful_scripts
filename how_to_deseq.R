library(DESeq2)
library(data.table)
library(tidyverse)
# First we are going to load in the functions that I've written as helper scripts
create_feature_path = "~/Documents/GitHub/helpful_scripts/create_feature_count_table.R"
make_deseq_path = "~/Documents/GitHub/helpful_scripts/make_deseq_dfs.R"
#you'll want to adjust the file paths accordingly
#source will bring the functions in these Rscripts into the current environment
source(create_feature_path)
source(make_deseq_path)

# this function creates a table wtih the first column being Geneid, then each
#column being the count in a sample, and the final column being the gene name
#it takes as an input the folder path where all the feature_counts files are
#the prefix is somethign which it appended by feature counts and the suffix
#is the bam suffix. Typically the bam suffix if it was aligned through
#our most recent version of the pipeline will be .Aligned.sorted.out.bam
#this table will be the input to our next function
murphy_feature = create_feature_count_table("/Users/annaleigh/Documents/data/deseq2_tutorial/murphy_feature_counts/",
                                          prefix = "X.SAN.vyplab.alb_projects.data.public_ribo_tags.murphy_royal_2020.STAR_aligned.",
                                          suffix = ".Aligned.sorted.out.bam")


#this next function 'make_deseq_df' returns a list of 2 formated data frames
#DESeq2 takes as input a count table and a table of metadata
#the 2 dataframes from this function are called 'conv_df' and 'coldata'
#conv_df because I've converted it into a format that DESeq2 likes
#and 'coldata' because it's metadata about the columns
make_deseq_dfs(murphy_feature)     #please note this does not assign to anything so the result will not be saved in the the environment
#I've made it easier to manipulate the data by providing
#additional arguments for what you want to grep (or serach) on the 
#column name as your 'baseline' and 'contrast' conditions
make_deseq_dfs(murphy_feature, base_grep = "IN_control", contrast_grep = "IP_control") #please note this does not assign to anything so the result will not be saved in the the environment
#if you just fill out the baseline, and there's only 2 it assumes that the 
#other is the contrast
make_deseq_dfs(murphy_feature, base_grep = "IN_control") #please note this does not assign to anything so the result will not be saved in the the environment
#you can assign the column and count matrixes to 2 separate variables using 
#using the dollar sign
murphy_counts = make_deseq_dfs(murphy_feature,base_grep = "IN_control", contrast_grep = "IP_control")$conv_df
#these functions print which columns they're keeping from the feature_count table to the console
#as a debugging
#you can select fewer with the command
#grep pattern
#for example if you only wanted samples 2,3,4
murphy_counts = make_deseq_dfs(murphy_feature,
                               grep_pattern = "2|3|4",
                               base_grep = "IN_control", 
                               contrast_grep = "IP_control")$conv_df

#but we don't want that for now
murphy_counts = make_deseq_dfs(murphy_feature,
                               base_grep = "IN_control", 
                               contrast_grep = "IP_control")$conv_df #note the $ we're only take one item of the list


murphy_meta = make_deseq_dfs(murphy_feature, 
                             base_grep = "IN_control", 
                             contrast_grep = "IP_control")$coldata #note the $ we're only take one item of the list


# When you do DESeq2, you should remove genes with very low counts,
#it can mess up the way DESeq2 normalizes counts and give you funky p-valeus
#here I wrote another helper function to do it
murphy_counts = filter_count_table(murphy_counts)

#This helper function will create another column using
#the 'cond' column of the metadata DF and give it whatever name you find
#more meaningful than 'base' and 'contrast'
#it also turns it into a factor with the 'base' condition as the first level
#the new column will be named 'comparison'
murphy_meta

murphy_meta = rename_relevel_for_deseq(murphy_meta,
                                       baseName = "input",
                                       contrastName = 'astrocyte')


murphy_meta
#now that we've done all the data sorting we actually get to the DESeq2 part
#first we make an object using the counts, the meta_data, and the column we
#want to compare. You should read the DESeq2 manual for
#more complicated design explanations
dds_murphy = DESeqDataSetFromMatrix(murphy_counts,
                                  colData = murphy_meta, 
                                  design = ~ comparison)
#that just created the object
#to actually run the analysis we just call "DESeq"
dds_murphy = DESeq(dds_murphy)

#we can quickly view the results with the function results()
results(dds_murphy)
#and then we can get a summary with summary on results
results(dds_murphy) %>% summary()
#we can pull the results data frame out and use the feature_counts table to 
#append the gene_names back on like this

res_murphy = results(dds_murphy) %>% 
    as.data.frame() %>% 
    rownames_to_column('Geneid') %>% 
    left_join(murphy_feature %>% dplyr::select(Geneid,gene_name))


