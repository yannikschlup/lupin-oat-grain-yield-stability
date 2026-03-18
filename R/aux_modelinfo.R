# Author: Lukas Graz             Date: 2024-10-10
# -----------------------------------------------
# given a Dgg & y: gives auxiliary infos about the model:
# y ~ treatment + ... + (1|Block_ID)
# - balancedness
# - anova of "..."
# - main & interaction effects of "..."
# -----------------------------------------------

aux_modelinfo <- function(D, hyp){
  stopifnot(nrow(D) == length(get(hyp$subset)))
  Dgg <- subset(D, get(hyp$subset) & !is.na(D[[hyp$varname]]))
  
  # first check if "year", "season", "location" have two levels in the data  
  use.year     <- if (length(unique(Dgg$year))     >= 2) TRUE else FALSE
  use.season   <- if (length(unique(Dgg$season))   >= 2) TRUE else FALSE
  use.location <- if (length(unique(Dgg$location)) >= 2) TRUE else FALSE

  # create y ~ 0 + treatment + treatment:`VARs` + `VARs-interactions` + (1|block_ID)
  # (remove the term if, use.year, use.season, use.location is FALSE)
  ff <- paste(hyp$varname, " ~ 0 + treatment", 
    if(use.year)     " + treatment:year",
    if(use.season)   " + treatment:season",
    if(use.location) " + treatment:location",
    if(use.year   & use.season)   " + year:season",
    if(use.year   & use.location) " + year:location",
    if(use.season & use.location) " + season:location",
    if(use.year   & use.season   & use.location) " + year:season:location",
    " + (1|block_ID)")

  # fit model
  fit <- lmer(as.formula(ff), data = Dgg)
  
  # COEFFICIENTS TABLE
  betas <- coef(summary(fit))[,"Estimate"]
  betas_t <- betas[grepl("treatment", names(betas))]
  betas_other <- betas[!grepl("treatment", names(betas))]
  trts <- grep("^treatment[0-9]+$", names(betas), value=TRUE) |> gsub("treatment", "", x=_)
  betas_t
  col_names <- c("reference", gsub("^treatment[0-9]+", "", names(betas_t)) |> unique() |> _[-1])
  COEFS <- matrix(c(
      betas_t[grepl("^treatment[0-9]+$", names(betas_t))],
      if (use.season)   betas_t[grepl("treatment.*season", names(betas_t))],
      if (use.location) betas_t[grepl("treatment.*location", names(betas_t))],
      if (use.year)     betas_t[grepl("treatment.*year", names(betas_t))]
    ), 
    dimnames = list(
      paste0("treatment", trts),  col_names), 
    ncol = length(col_names)
  ) |> zapsmall()
  COEFS <- rbind(
    `mean(collumn)` = apply(COEFS, 2, mean),
    `sd(collumn)` = apply(COEFS, 2, sd),
    COEFS
  )
  # @lukas add 3way interactions

  COEFS

  aux_modelinfo <- list(
    COEFS = COEFS,
    betas = betas,
    betas_other = betas_other,
    anova = anova(fit),
    fit = fit
  )
  class(aux_modelinfo) <- "aux_modelinfo"
  aux_modelinfo
}

aux_modelinfo_cache <- function(D, hyp){
  xfun::cache_rds({
      aux_modelinfo(D, hyp)
    },
    hash = list(
      hyp,
      D,
      as.character(functionBody(aux_modelinfo)),
      formals(aux_modelinfo)
    ),
    dir = "cache/",
    file = paste0("aux_modelinfo_", 
                  hyp$H_ID, "_"
                  )
  )
}


print.aux_modelinfo <- function(x, ...){
  is.knitr <- !is.null(getOption("knitr.in.progress"))
  if(is.knitr) cat("\n\n```{r}\n")

  cat("\n\nAnova: \n")
  print(x$anova)

  cat("\n\nCoefficients from model:  \ny ~ 0 + trt + trt:`VARs` + `VARs-interactions` + (1|block_ID):  \n")
  print((x$COEFS))
  cat("\n\nRemaining (non-treatment) Coefficients: \n")
  print(x$betas_other)

  if(is.knitr) cat("\n```\n")
  invisible(x)
}



if(FALSE){
  aux_modelinfo(subset(D, get(hyp$subset)), hyp)
}




  # attr(..., "NA_unbalanced") <-   tryCatch(xtabs(~cv_oat + cv_lupin, data = Dgg, addNA = TRUE))

    # cat("\n\nWhich Plots were NAs (removed): \n")
  # print(paste(attr(p, "whichNAs"), sep=", "))
  # if (length(attr(p, "NA_unbalanced")) > 5){
  #   cat("\n\nxtable of NAs: \n")
  #   print(knitr::kable(attr(p, "NA_unbalanced")))
  # }
