
# human = download_appris('human')

download_appris = function(species){
    #read in the appris annotations on principal transcripts
    if(species == "human"){
        appr_anno= fread("http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/appris_data.principal.txt",header = F)
    }else if(species == "mouse"){
        appr_anno = fread("http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/mus_musculus/GRCm38/appris_data.principal.txt",header = F)
    }else{
        print("Species currently either 'human' or 'mouse'")
    }
    # add a column that looks for if the appris transcript is called a 'principal' 
    # or a alternative

    appr_anno[grep("PRINCIPAL",V5), Principal := as.integer(str_extract(V5,'[[:digit:]]+'))]
    appr_anno[!is.na(Principal), NumPrincipal := .N, by = V2]
    return(appr_anno)
}

show_appris_transcripts <- function(gene, ident = "hugo_symbol",species = "human",appr_anno = ""){
    if (!require("pacman")) install.packages("pacman")
    library(pacman)
    p_load("stringr","data.table")
    # it takes some time to download the appris annotation so we can load it into the database as the start if we like
    # otherwise we download the appropiate one

    suppressWarnings(if(appr_anno == ""){
        appr_anno = download_appris(species)
    })

    #if the identifer is hugo symbol, check to see if you can find it
    if(ident == "hugo_symbol"){
        if(!gene %in% appr_anno[,V1]){
            return("Hugo Symbol not found, try an alias or searching Ensembl gene")
        } else{
            return(appr_anno[V1 == gene])
        }
    } else if(ident == "ens"){
        if(!gene %in% appr_anno[,V2]){
            return("Ensembl gene not found")
        } else{
            return(appr_anno[V2 == gene])
        }
    }
}

choose_principal_transcript = function(gene, ident = "hugo_symbol",species = "human",appr_anno = ""){

    appris_subset = show_appris_transcripts(gene, ident = ident,species = species,appr_anno)
    # okay the first case is that a gene only has one row with appris "Principal" not NA
    if(nrow(appris_subset[!is.na(Principal)]) == 1){
        principal = appris_subset[!is.na(Principal),V3]
    }else(
        browser()
    )
    return(principal)
}
