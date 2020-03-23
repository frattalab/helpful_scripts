# Basic function to  human to mouse gene names and vice versa
#give genes as a list , from species is what gene ids you have and to species is what you're going to
convertGeneListSpecies <- function(genelist, from_species = "hsapiens_gene_ensembl", to_species = "mmusculus_gene_ensembl"){
  if (!require("pacman")) install.packages("pacman")
  library(pacman)
  p_load(biomaRt)
  #listDatasets(mart = useMart("ensembl")) if you run this line it will  give you all the datasets you can choose from
  ensembl.from = useMart("ensembl", dataset = from_species)
  ensembl.to = useMart("ensembl", dataset = to_species)
  mapped_list = getLDS(attributes=c("ensembl_gene_id"),
                       filters="ensembl_gene_id", values=genelist, mart=ensembl.from,
                       attributesL=c("ensembl_gene_id"), martL=ensembl.to)

  return(mapped_list)
}
#function coverts two different types of genes; using the annotation dbi
convertGeneListTypes = function(genelist, from_type = "ENTREZID", to_type = "ENSEMBL", species = "mmusculus",print_keytypes = "FALSE"){
  #ensure you have
  if (!require("pacman")) install.packages("pacman")
  library(pacman)
  p_load(AnnotationDbi)
  # if we just wanted to prrint key type to check
  if(print_keytypes == TRUE){
    print(keytypes(org.Mm.eg.db))
    return()
  }
  # if we tried to type in a species that doesn't exists because I'm a bad speller
  if(!species %in% c("mmusculus","hsapiens")){
    print("Species must be 'mmusculus' or 'hsapiens'")
  }
  
  
  
  if( species == "mmusculus"){
    
      anno_db = p_load(org.Mm.eg.db)
      mapped_list = as.data.table(AnnotationDbi::select(org.Mm.eg.db, keys =  unique(genelist), keytype = from_type, columns = c("SYMBOL",to_type)))
  
  }else if(species == "hsapiens"){
    anno_db = p_load(org.Hs.eg.db)
    mapped_list = as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys =  unique(genelist), keytype = from_type, columns = c("SYMBOL",to_type)))
    }

  return(mapped_list)

}
