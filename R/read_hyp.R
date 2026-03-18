# Author: Lukas Graz             Date: 2024-10-14
# -----------------------------------------------
# Reads hypotheses from hypotheses.csv and does all nice checks
# -----------------------------------------------

# source("R/DataPrep.R")
hypotheses <- read.csv("./hypotheses.csv", strip.white = TRUE, skip = 2)

# make sure if data loaded
if(!exists("D_typed"))
  source("R/DataPrep.R")
stopifnot(exists("D_typed"))

# check if all varnames are in the data
stopifnot(all(hypotheses$varname %in% colnames(D)))

# check if value in subset is a variable in ls()
stopifnot(all(hypotheses$subset %in% ls(), na.rm = TRUE))

# check if H_ID and varname are all non-na's
stopifnot(all(!is.na(hypotheses$H_ID)))
stopifnot(all(!is.na(hypotheses$varname)))

# check if H_ID and varname are always (nonempty) strings
stopifnot(all(nchar(hypotheses$H_ID) > 0))
stopifnot(all(nchar(hypotheses$varname) > 0))

# check if H_ID is unique
stopifnot(length(unique(hypotheses$H_ID)) == nrow(hypotheses))



# replace NA by default values
replace_na_in_target_by_source <- function(target, source, value=NULL){
  is.na.ind <- is.na(hypotheses[[target]])
  if (any(is.na.ind)) 
    hypotheses[is.na.ind, target] <<- if(is.null(value)) hypotheses[is.na.ind, source] else value
}
replace_na_in_target_by_source("ggyLabel", "varname")
replace_na_in_target_by_source("plotTruncate", value=FALSE)
replace_na_in_target_by_source("robustify", value=FALSE)
replace_na_in_target_by_source("transform", value="identity")

rownames(hypotheses) <- hypotheses$H_ID



# -----------------------------------------------
# defines transformations for data
# -----------------------------------------------

library(DescTools)

trafo      <- function(x, transform = "identity"){
  if(is.na(transform) | transform == "identity") x else 
    if(transform == "log10") {stopifnot(all(x>0, na.rm=TRUE));  log10(x)} else 
    if(transform == "logst") {stopifnot(all(x>0, na.rm=TRUE));  LogSt(x)} else 
    if(transform == "sqrt")  {stopifnot(all(x>=0, na.rm=TRUE)); sqrt (x)} else 
      stop("transform not implemented")
}

inv_trafo <- function(x, transform = "identity"){
  if(is.na(transform) | transform == "identity") x else 
    if(transform == "log10") 10^x else 
    if(transform == "logst") LogStInv(x) else 
    if(transform == "sqrt" ) {stopifnot(all(x>=0, na.rm=TRUE)); x^2} else 
      stop("transform not implemented")
}
for (f in c("log10", "logst", "sqrt")) 
  stopifnot(cor(1:10, inv_trafo(trafo(1:10, f), f)) > 0.999)



# transform varname according to transform (for unique variables)
stopifnot(exists("D_typed"))
if (exists(".alrdy_transformed")) {
  .alrdy_transformed <- TRUE
  message("Data already transformed, skipping...")
} else {
  duplicated_varnames <- duplicated(hypotheses$varname)
  for (i in seq_along(hypotheses$varname)){
    if (hypotheses$transform[i] %in% c("identity", NA)){
      next
    } else if (hypotheses$transform[i] %in% c("log10", "logst", "sqrt")){
      old_varname <- hypotheses$varname[i]
      new_varname <- paste0(hypotheses$transform[i], "_", hypotheses$varname[i])
      hypotheses$varname[i] <- new_varname
      if (!duplicated_varnames[i]) 
        D[[new_varname]] <- trafo(D[[old_varname]], hypotheses$transform[i])
    } else {
      stop("transform not implemented")
    }
  }
}



