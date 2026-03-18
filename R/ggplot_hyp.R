# Author: Lukas Graz             Date: 2024-10-07
# -----------------------------------------------
# infrastructure to plot data by hypothesis
# -----------------------------------------------


# @lukas: argue in methods to not test contrast (nor visualize p-values in plots) to mitigate multiple testing
# - [x] lukas: "Pure lupins:"
# - [ ] lukas: Lendend: mixed/pure (only if not only mixed)
# - [x] lukas: include outliers but color red

# @END: set all titles / captions /...  in the correct place

library(ggplot2)
library(lmerTest, quietly = TRUE)
if(!require(remef)) {
  devtools::install_github("hohenstein/remef")
  library(remef)
}

theme_set(
  theme_gray(base_family = "Times")
)

# custom labeller
custom_labeller <- function(value) {
  ifelse(value == "pure", "Pure lupins: ", paste(value, "mixed with:", sep=" "))
}

# truncate outliers
truncate_y <- function(y){
  iqr <- IQR(y, na.rm = TRUE)
  cencor_IQR <- quantile(y, c(0.25, 0.75), na.rm = TRUE) + c(-1, 1) * if(iqr != 0) iqr * 1.5 else Inf
  # cencor_5pers <- quantile(y, c(0.05, 0.95), na.rm = TRUE)
  cencor_5pers <- c(Inf, -Inf)
  lwr <- min(cencor_IQR[1], cencor_5pers[1], na.rm = TRUE)
  upr <- max(cencor_IQR[2], cencor_5pers[2], na.rm = TRUE)
  truncated <- y < lwr | y > upr
  truncated <- factor(truncated, levels = c(FALSE, TRUE), labels = c("no", "yes"))
  y_truncated <- y
  if (any(y < lwr, na.rm = TRUE))
    y_truncated[y < lwr] <- lwr
  if (any(y > upr, na.rm = TRUE))
    y_truncated[y > upr] <- upr

  list(
    y_truncated = y_truncated,
    truncated = truncated,
    lwr = lwr,
    upr = upr
  )
}
# truncate_y(c(1:20, NA, 10000))

ggplot_hyp <- function(hyp, D){
  set.seed(123)

  # DATA PREP =========

  # Subset data
  Dgg <- subset(D, get(hyp$subset))
  Dgg$y <- Dgg[[hyp$varname]]

  # remove NA's
    NAs <- which(is.na(Dgg$y))
    if (length(NAs) > 0){
      message(hyp$varname, ": remove H_ID(s) with NA's: ", paste(as.character(Dgg$H_ID[NAs]), collapse=", "))
      Dgg <- Dgg[-NAs,]
    }

  # Fix levels
    # 1. oat
    levels(Dgg$cv_oat) <- c(levels(Dgg$cv_oat), "pure")
    Dgg$cv_oat[is.na(Dgg$cv_oat)] <- "pure"
    Dgg$cv_oat <- droplevels(Dgg$cv_oat)
    # 2. lupin
    levels(Dgg$cv_lupin) <- c(levels(Dgg$cv_lupin), "nothing")
    Dgg$cv_lupin[is.na(Dgg$cv_lupin)] <- "nothing"
    Dgg$cv_lupin <- droplevels(Dgg$cv_lupin)

  # Truncation
    Dgg$truncated <- truncate_y(Dgg$y)$truncated
    # truncate outliers to max/min of cencor
    if (hyp$plotTruncate)
      Dgg$y <- truncate_y(Dgg$y)$y_truncated  # TODO: uncomment me

  # # Remove Blockeffect
  #   if (hyp$plotAdj){
  #     trt_fit <- lmer(y ~ 0 + treatment + (1|year) + (1|season) + (1|location) + (1|block_ID), data = Dgg)
  #     # keep only the treatment effects:
  #     Dgg$y <- keepef(trt_fit, fix = rownames(coef(summary(trt_fit))), keep.intercept = FALSE)
  #   }
  # DATA PREP END =================


  # PLOT =========
  mixed_and_pure_present <- length(levels(droplevels(Dgg$stand)))> 1
  truncated_present <- hyp$plotTruncate && length(levels(droplevels(Dgg$truncated)))> 1
  # basic plot
  p <- ggplot(Dgg, aes_string(
    y = "y",
    x = "cv_lupin",
    fill = "cv_lupin",
    linetype = if (mixed_and_pure_present) "stand"
    )) +
    labs(title = hyp$ggTitle, 
        y = hyp$ggyLabel, 
        tag=if(!is.na(hyp$ggTag))hyp$ggTag, 
        x = NULL, 
        caption = if(!is.na(hyp$ggCaption))hyp$ggCaption, 
        subtitle = hyp$ggSubtitle, 
        alt=hyp$H_ID, 
        label=hyp$H_ID)+ 
    # scale_y_continuous(transform = hyp$transform) +
    facet_wrap(~cv_oat, nrow=1, strip.position = "bottom", labeller = as_labeller(custom_labeller))

  # add plot_type specific layers
  plot_type <- hyp$plot_type
  p <- if (plot_type == "box") {
    p + geom_boxplot(outliers =FALSE, notch = TRUE) +
      geom_jitter(height=0, alpha=0.3, if(truncated_present) aes(col = truncated))
  } else if (plot_type == "likert"){
    p + geom_violin() +
      geom_jitter(height=0, alpha=0.3, if(truncated_present) aes(col = truncated))
  } else {
    stop("plot_type not implemented")
  }

  # Legend:
  p <- p + theme(legend.position = "bottom") +
    guides(fill = "none") +
    if(truncated_present) scale_color_manual(values = c("no" = "black", "yes" = "#dd1b62"))
    
    if (mixed_and_pure_present) p <- p + scale_linetype_manual(values = c("mixed" = "solid", "pure" = "twodash")) # Solid for mixed, twodash for pure

  class(p) <- c("ggplot_hyp", class(p))
  attr(p, "hyp") <- hyp
  p
}

# TESTS
if (FALSE){
  source("R/DataPrep.R")
  # source("R/ggplot_hyp.R")
  source("R/Models.R")
  source("R/aux_modelinfo.R")
  source("R/read_hyp.R")
  source("R/eval_hyp.R")

  
  hyp <- hypotheses[8,]
  hyp <- hypotheses[16,]
  hyp <- hypotheses["H9",]
  hyp <- hypotheses["H3.l",]
  hyp <- hypotheses["H8.lo",]
  hyp <- hypotheses["H1.ol.fall",]
  hyp <- hypotheses["H1.ol.spring",]
  hyp$transform <- "identity"
  hyp$transform <- "atanh"
  hyp$transform <- "log10"
  ggplot_hyp(hyp, D)
}


# CACHE wrapper
ggplot_hyp_cache <- function(hyp, D){
  xfun::cache_rds({
      ggplot_hyp(hyp, D)
    },
    hash = list(
      hyp,
      D,
      as.character(functionBody(ggplot_hyp)),
      formals(ggplot_hyp)
    ),
    dir = "cache/",
    file = paste0("ggplot_hyp_", 
                  hyp$H_ID, "_"
                  )
  )
}

#' Removing effects I have reinvented the wheel
#' use package `remef` instead
#' the following illustrates that we get the same results
if (FALSE){
  # setup
  hypotheses <- read.csv("./hypotheses.csv")
  hyp <- hypotheses[9,]
  plotAdj = TRUE

  # old function code:
  {
    set.seed(123)

    # Subset data
    Dgg <- subset(D, get(hyp$subset))
    Dgg$y <- Dgg[[hyp$varname]]

    # remove NA's
    NAs <- which(is.na(Dgg$y))
    if (length(NAs) > 0){
      message(hyp$varname, ": remove plot_ID(s) with NA's: ", paste(as.character(Dgg$plot_ID[NAs]), collapse=", "))
      Dgg <- Dgg[-NAs,]
    }

    # remove effect from y of year, season, location, and block using mixed model
    if (plotAdj){
      trt_fit <- lmer(y ~ 0 + treatment + (1|year) + (1|season) + (1|location) + (1|block_ID), data = Dgg)
      Trt_signal <- coef(summary(trt_fit))[,"Estimate", drop=FALSE]
      trt_nr <- as.numeric(gsub("treatment", "", rownames(Trt_signal)))
      TrtEstimate <- c()
      TrtEstimate[trt_nr] <- Trt_signal[,"Estimate"]
      Dgg$y_clean <- residuals(trt_fit) + TrtEstimate[Dgg$treatment]
    }
  }

  # comparison with remef 
  library(remef)
  y_keepef <- keepef(trt_fit, fix = rownames(coef(summary(trt_fit))), keep.intercept = FALSE)
  pairs(data.frame(
    y = Dgg$y,
    y_clean = Dgg$y_clean,
    keepef = y_keepef
  ))
  cor(Dgg$y_clean, y_keepef)
}
