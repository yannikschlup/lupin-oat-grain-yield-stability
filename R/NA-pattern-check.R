# Author: Lukas Graz             Date: 2024-10-02
# -----------------------------------------------
# Check if the are false non-NA's
# -----------------------------------------------


### NA Pattern:
# Using the `Applies to` in the metadata, we to expect NAs 
# `bad_nonNAs_positions` gives us all positions where we expect values but got NAs
get_bad_nonNAs_positions <- function(v, applies_to){
  if (applies_to == "all"){ 
    v[] <- NA
  }else if (applies_to %in% c("mixed", "mixture")) {
    v[is.mixed] <- NA
  }else if (applies_to %in% c("oat")){
    v[is.oat] <- NA
  }else if (applies_to %in% c("lupin")){
    v[is.lupin] <- NA
  }else 
    stop("applies_to not recognized")
  which(!is.na(v))
}

bad_nonNAs_positions <- list()
stopifnot(which(is.na(D_meta$`Applies to`)) == 1:11)
for (clmn in names(D_typed)[-(1:11)]){
  bad_nonNAs_positions[[clmn]] <- get_bad_nonNAs_positions(D_typed[[clmn]], D_meta[clmn==names(D_typed),"Applies to"])
}
lapply(bad_nonNAs_positions, length) |> unlist() |> cbind()
bad_nonNAs_positions$oat_tkw # corrupt index
D_typed$oat_tkw[464] # value
