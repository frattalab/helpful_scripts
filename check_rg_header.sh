#Script to check if @RG flag is missing at start of header line in BAM files and edit accordingly 

#Some older lab data is missing @RG flag at the start of line
#Missing flag was accepted by older versions of BAM file tools (e.g. samtools, pysam etc.) - newer versions now refuse to read/open these files
#Adding flag makes BAM file compatible with latest versions of these BAM file tools

## REQUIREMENTS
#samtools version 1.3.1 (slightly more recent versions may well work but I haven't tried them (problem came up with samtools 1.9)

## STEPS
#1. Find all BAM files in current directory (and any subdirectories contained withing working dir)
#2. Check if RG header line begins with @RG
#3. If @RG missing then add to header line to create reformatted BAM file
#4. Index reformatted BAM files
#5. (AFTER SCRIPT HAS RUN) check reformatted files look ok then replace old files with reformatted files (mv <reformatted>.bam <old_file>.bam)
#(could probably do step 5 within this script but I'm a bit paranoid about losing the original data...)

#NOTE
#Script works as is because there is only 1 read group (and so 1 header) in the BAM file (for me M323K & F210I BAMs).
#If BAM files processed with RNA_seq_pipeline_v8.sh then this is probably the case but is worth double checking before running the script

#Activate conda environment with samtools 1.3.1 
source activate samtools_old_test

#find BAM files in working dir

find . -name '*.bam' > bam_file_list.txt

#cat bam_file_list.txt

line_no=$(wc -l < bam_file_list.txt)

#loop through each identified bam file in wd and subdirectories 

for x in $(seq 1 $line_no); do
    BAM=$(awk -v val="$x" '{if(NR==val) print $0}' bam_file_list.txt)
        
    #check if RG header line begins with @RG
    #Print filepaths of ok & to fix files to separate text files
    samtools view -h $BAM | head | awk '{if(FNR==2) print $1}' | awk -v name="$BAM" '{if($1=="@RG") {print name >> "bam_rg_ok_list.txt";next} else {print name >> "bam_rg_to_fix_list.txt"}}'

done

#Add @RG to header line of files in to_fix text file 

line_no=$(wc -l < bam_rg_to_fix_list.txt)

for x in $(seq 1 $line_no); do 
    input_BAM_path=$(awk -v val="$x" '{if(NR==val) print $0 }' bam_rg_to_fix_list.txt)    
    BAM_name=$(awk -v val="$x" '{if(NR==val) print $0 }' bam_rg_to_fix_list.txt | awk -F. '{ OFS = "." } NF{NF-=1} {print $0}')

    samtools view -h $input_BAM_path | awk 'NR==2 { print "@RG\t" $0 } NR!=2 {print}' | samtools reheader - $input_BAM_path > ${BAM_name}_rg_fixed.bam
    samtools index -b ${BAM_name}_rg_fixed.bam > ${BAM_name}_rg_fixed.bam.bai
done
