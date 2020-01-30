
get_principal_transcript <- function(gene, ident = "hugo_symbol",species = "human"){
    if (!require("pacman")) install.packages("pacman")
    library(pacman)
    p_load("stringr","data.table")
    #read in the appris annotations on principal transcripts
    if(species == "human"){
        appr_anno= fread("http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/appris_data.principal.txt",header = F)
    }else if(species == "mouse"){
        appr_anno = fread("http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/mus_musculus/GRCm38/appris_data.principal.txt")
    }else{
        print("Species currently either 'human' or 'mouse'")
    }
    
    appr_anno[grep("PRINCIPAL",V5), Princ := as.integer(str_extract(V5,'[[:digit:]]+'))]

    #if the identifer is hugo symbol, check to see if you can find it
    if(ident == "hugo_symbol"){
        if(!gene %in% appr_anno[,V1]){
            return("Hugo Symbol not found, try an alias or searching Ensembl gene")
        } else{
            print(appr_anno[V1 == gene])
        }
    } else if(ident == "ens"){
        if(!gene %in% appr_anno[,V2]){
            return("Ensembl gene not found")
        } else{
            print(appr_anno[V2 == gene])
        }
    }
}
