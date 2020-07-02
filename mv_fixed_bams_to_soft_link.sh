#!/usr/bin/env bash

##Bash script to move BAM  to their soft linked location
#should be ls <command> where command would select the files of choice
#e.g. searching for all files in current directory, would pass *rg_fixed* IN QUOTES
rg_ls_command=$1

#echo $rg_ls_command
#if file has content it empties it - if not creates empty file
#don't want to add commands over and over
truncate -s 0 fixed_bams_to_move_list.txt


 ls ${rg_ls_command} | for i in $(cat -); do
   fixed_file_name=$(echo $i | cut -d. -f1) #get rid of exrention
   original_name=${fixed_file_name%"_rg_fixed"} # strip _rg_fixed suffix from bam name - gets name of original file
   ext=$(echo $i | awk -F "." '{for (i=2; i<=NF; i++) {if ($i !="") {printf "."$i};}}') #gets extension of file in question (up to two e.g. .bam.bai)
   #https://stackoverflow.com/questions/8984720/awk-joining-n-fields-with-delimiter
   #echo 'fixed file name is' ${fixed_file_name}${ext}

   #original name should be a soft link to the original location - want to move fixed file to this location
   #if ll <original_name> , the last field returned provides path to which soft link points
   #.bam.bai for original files are not linked - need to add .bai to path to corresponding .bam file
   if [[ "$ext" == ".bam" ]]
   then

     original_location=$(ls -l ${original_name}${ext} | awk '{print $NF}')
     #echo 'original location is' ${original_location}
     echo ${i} ${original_location} >> fixed_bams_to_move_list.txt

   elif [[ $ext == ".bam.bai" ]]
   then

     original_bam_location=$(ls -l ${original_name}.bam | awk '{print $NF}')
     original_index_location="${original_bam_location}.bai"
     #echo 'original index location is' $original_index_location
     echo ${i} ${original_index_location} >> fixed_bams_to_move_list.txt
   fi

 done

echo 'the following commands will be performed'
#https://stackoverflow.com/questions/7551991/add-a-new-column-to-the-file
awk 'BEGIN{FS=OFS=" "}{print "mv" OFS $0}' fixed_bams_to_move_list.txt

#function to loop through fixed_bams_to_move_list - mv element in first column to element in second
function mv_files {
  echo 'file containing files to move is' $1
  while IFS=" " read -r line; do
    #printf "line in text from file is:\n${line}\n"
    mv ${line}
    printf "\nFollowing command has been completed:\nmv ${line}\n"
  done < $1
}

#user prompt to confirm whether to perform all moves
#https://stackoverflow.com/questions/226703/how-do-i-prompt-for-yes-no-cancel-input-in-a-linux-shell-script/27875395#27875395
echo 'Do you wish to proceed with all the above commands (put 1 or 2 to select yes or no)?'
select yn in "yes" "no"; do
  case $yn in
    yes ) echo 'you selected to perform these commands'; mv_files fixed_bams_to_move_list.txt;
     break;;
    no ) echo 'You have opted to not perform the above commands. Aborting...'; exit;;
  esac
done
