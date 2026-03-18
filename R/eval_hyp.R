# Author: Lukas Graz             Date: 2024-10-10
# -----------------------------------------------
# evaluates (= model + plot) a hypothethis (with caching)
# -----------------------------------------------

source("R/ggplot_hyp.R")
source("R/Models.R")
source("R/aux_modelinfo.R")

eval_hyp <- function(hyp, D){
  # @lukas: check if hyp is valid
  L <- list(
    hyp = hyp,
    model = eval_model(hyp, D),
    plot = if (hyp$plot_type != "") ggplot_hyp_cache(hyp, D)
    )
  class(L) <- "eval_hyp"
  L
}

print.eval_hyp <- function(x){
  is.knitr <- !is.null(getOption("knitr.in.progress"))
  cat("\n\n## ",x$hyp$H_ID,": ", x$hyp$varname, "\n\n")
  # PLOTS
  if (!is.null(x$plot)){
    cat("\n\n**Plots**\n\n")
      suppressMessages(
        print(x$plot)
        )
    # cat("\n\n**Aux Model Info**\n\n")
    # print(x$aux_modelinfo)
  }

  # MODEL
  if (!is.null(x$model)){
    cat("\n\n**Model: ", x$hyp$model_type ,"**\n\n")
    print(x$model)

    cat("\n\n**Model Summary:**\n\n")
    if(is.knitr) cat("\n\n```text\n")
    suppressMessages(print(summary(x$model$model)))
    if(x$hyp$robustify){
      cat("Robustness weights from robustlmm::rlmer\n")
      summarizeRobWeights(x$model$weights$wgt.e, digits = 3, header = "\nRobustness weights for the residuals:")
      summarizeRobWeights(x$model$weights$wgt.b, digits = 3, header = "\nRobustness weights for the random effects:")
    }
    if(is.knitr) cat("\n```\n")
  }

  # HYP
  cat("\n\n**Hypothesis:**\n\n")

  # if(is.knitr) cat("\n\n```markdown\n")
  # print(t(x$hyp))
  
  # remove empty rows
  Thyp <- t(x$hyp)
  rm_ind <- is.na(Thyp) | Thyp == "" | Thyp == "NA" | row.names(Thyp) == "H_ID"
  print(knitr::kable(Thyp[!rm_ind,,drop=FALSE], format = "markdown"))
  # if(is.knitr) cat("\n```\n")
}


# TESTS
if (FALSE){
  source("R/DataPrep.R")
  source("R/read_hyp.R")
  i=5
  for (i in 1:nrow(hypotheses)){
    print(i)
    Sys.sleep(1)
    hyp <- hypotheses[i,]
    x <- eval_hyp(hypotheses[i,], D) 
    print(x)
    # eval_model(hypotheses[i,], D)

  }
  eval_model(hypotheses["H5",], D)
  eval_model(hypotheses["H81",], D)
}
