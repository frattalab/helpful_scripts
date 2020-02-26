# -----takes a data.frame and either converts a specific column to a numeric
# ------- or find all the non numeric columns and converts them to numerics
# ----use liek this convert_to_numeric(klim_meta, "cell")
# ---- or like this convert_to_numeric(klim_meta)
convert_to_numeric <- function(dataframe){
        # ---select all the columns which are characters
        chararcter_columns = dataframe %>% 
            select(which(sapply(.,is.character))) %>% colnames()
        for (c in chararcter_columns){
            new_col = paste0(c,"_numeric")
            dataframe[, eval(quote(new_col))] = as.numeric(as.factor(dataframe[, eval(quote(c))]))
        }
return(dataframe)
}
      


