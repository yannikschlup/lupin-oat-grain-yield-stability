# Author: Lukas Graz             Date: 2024-10-07
# -----------------------------------------------
# List of models -> models are applied according to hypotheses.csv
# -----------------------------------------------

#' each model is a function that takes the arguments:
#' - subset
#' - data (per default D from R/DataPrep.R)
#' Returns a list with:
#' - model: the model object
#' - results,
#' - p-value(s) (non-adjusted)

library(lmerTest, quietly = TRUE)
library(robustbase)
library(DescTools)
library(robustlmm)
library(simr)
library(emmeans) #install.packages("emmeans")
source("R/rlmerTest.R")


# nicer names for contrasts
# c.f. https://stackoverflow.com/questions/24515892/r-how-to-contrast-code-factors-and-retain-meaningful-labels-in-output-summary
named.contr.sum <- function(x, ...) {
    if (is.factor(x)) {
        x <- levels(x)
    } else if (is.numeric(x) & length(x)==1L) {
        stop("cannot create names with integer value. Pass factor levels")
    }
    x<-contr.sum(x, ...)
    colnames(x) <- apply(x,2,function(x) 
         paste(names(x[x>0]), sep="-")
    )
    x
}

# options(contrasts = c("contr.sum", "contr.poly"))
options(contrasts = c("named.contr.sum", "contr.poly"))
message("All models are using 'contr.sum' contrasts (with nicer names)")

#' drops unavailable terms from " + year * season * location + (1|block_ID)"
get_str_formula_rest <- function(Dgg, hyp){
  # first check if "year", "season", "location" have two levels in the data  
  use.y <- if (length(unique(Dgg$y)) >= 2) TRUE else FALSE
  use.l <- if (length(unique(Dgg$l)) >= 2) TRUE else FALSE

  paste(
    # hyp$varname, 
    # " ~ 0 + trt", 
    # if(use.y) " + trt:y",
    # if(use.l) " + trt:l",
    if(use.y) " + y",
    if(use.l) " + l",
    if(use.y & use.l)  " + y:l",
    " + (1|block_ID)")
}

get_last_column_from_df <- function(results){
  ind <- !(rownames(results) %in% c("y", "l", "y:l","y2022", "lRe", "y2022:lRe")) 
  x <- results[[names(results)[ncol(results)]]]
  # nams <- rownames(x)
  # x <- drop(x)
  # names(x) <- nams
  names(x) <- rownames(results)
  x[ind]
}

get_weights_via_rlmer <- function(form, Dgg){
  rfit_tmp <- rlmer(form, data=Dgg)
  RANEF <- cbind(
    ranef(rfit_tmp)[[1]],
    wgt.b = robustlmm:::wgt.b(rfit_tmp),
    wgt.ranef = unname(ranef(rfit_tmp)[[1]] * robustlmm:::wgt.b(rfit_tmp)), 
    block_ID = rownames(ranef(rfit_tmp)[[1]])
  )
  # wgt.e * wgt.b
  list(
    wgt = robustlmm:::wgt.e(rfit_tmp) * merge(rfit_tmp@frame, RANEF, by="block_ID")[,"wgt.b"],
    wgt.e = robustlmm:::wgt.e(rfit_tmp),
    wgt.b = robustlmm:::wgt.b(rfit_tmp)
  )
}

.aux <- function(RHS_formula_core,
                 diagnostics_fun = NULL,
                 analysis_fun   = lmerTest:::anova.lmerModLmerTest,
                 get_p_vals_fun = \(df) try(get_last_column_from_df(df)),
                 formula_rest   = NULL) {

  function(hyp, D = D) {

    varname <- hyp$varname
    Dgg     <- subset(D, get(hyp$subset) & !is.na(D[[varname]]))

    if (is.null(formula_rest))
      formula_rest <- get_str_formula_rest(Dgg, hyp)

    form <- as.formula(paste0(varname, " ~ ", RHS_formula_core, formula_rest))

    #robust weights
    weights <- NULL
    if (hyp$robustify) {
      wgt     <- get_weights_via_rlmer(form, Dgg)
      weights <- wgt$wgt
    }

    #model fit, diagnostics, p‑values
    model       <- lmer(form, data = Dgg, weights = weights)
    diagnostics <- if (!is.null(diagnostics_fun)) diagnostics_fun(model)
    results     <- analysis_fun(model)

    e      <- tryCatch(get_p_vals_fun(results), error = identity)
    p_vals <- if (!inherits(e, "error")) e else NA

    emmeans_posthoc <- NULL

    if (hyp$model_type == "M3") {

      #cv contrasts
      if (hyp$subset == "is.oat") {
        emmeans_posthoc <- list(
          cultivar = pairs(emmeans(model, specs = ~cv_oat))
        )

      } else if (hyp$subset == "is.lupin") {
        emmeans_posthoc <- list(
          cultivar = pairs(emmeans(model, specs = ~cv_lupin))
        )

      } else if (hyp$subset == "is.mixed") {
        emmeans_posthoc <- list(
          cultivar = pairs(emmeans(model, specs = ~cv_oat:cv_lupin))
        )
      }

      #location contrasts, 2022 & 2023 (all five sites)
      if (all(c("2022", "2023") %in% model@frame$y)) {
        loc_22_23 <- pairs(
          emmeans(model,
                  specs = ~l,
                  by    = "y",
                  at    = list(y = c("2022", "2023")))
        )
        emmeans_posthoc$loc_22_23 <- loc_22_23
      }

      #location contrasts, 2024 (Gross / Novi / Vilbel only)
      if ("2024" %in% model@frame$y) {
        present_lvls <- intersect(c("Gross", "Novi", "Vilbel"),
                                  levels(model@frame$l))
        if (length(present_lvls) >= 2) {
          loc_24 <- pairs(
            emmeans(model,
                    specs  = ~l,
                    at     = list(y = "2024"),
                    params = list(levels = list(l = present_lvls)))
          )
          emmeans_posthoc$loc_24 <- loc_24
        }
      }
    }

    RETURN <- list(
      varname         = varname,
      model           = model,
      formula         = formula(model),
      diagnostics     = diagnostics,
      results         = results,
      weights         = if (hyp$robustify) wgt,
      p_vals          = p_vals,
      emmeans_posthoc = emmeans_posthoc
    )
    class(RETURN) <- "model.analysis"
    RETURN
  }
}


# given a fit: tests if `which_coef` is different than `val` (+Bonferroni)
pval_FUN_which_coef_diff_from_val <- function(val = 0){
  function(model){
    p <- length(fixef(model))
    Lmat <- diag(p)
    tab <- lmerTest:::rbindall(lapply(1:p, function(i) contest1D(model, 
        Lmat[i, ], rhs = val)))
    rownames(tab) <- names(fixef(model))
    colnames(tab)[4] <- paste0("t=(Est-",val,")/Std.Err")
    tab[1:9,]
    # COEFs <- COEFs[grep(grep_include, rownames(COEFs)),]
    # estimates <- COEFs[, "Estimate"]
    # which_row <- names(which_coef(estimates - val))
    # COEF <- COEFs[which_row, , drop=TRUE]
    # t.value <- (COEF[["Estimate"]] - val) / COEF[["Std. Error"]]
    # p.value <- 2 * pt(-abs(t.value), df = COEF[["df"]]) # two-sided
    # # bonferroni correction:
    # p.value <- p.adjust(p.value, method = "bonferroni", n = nrow(COEFs))
    # names(p.value) <- paste0("adj. p-val H: |",which_row," - ", val, "| == 0")
    # p.value
  }
}

#' create contrast vector for A - B  (weighted) givn A_levels and B_levels
contrast_vector_A_minus_B <- function(fit, A_levels, B_levels, trt_name = "trt"){
  diff_contrast <- vector("numeric", length = length(fixef(fit)))
  names(diff_contrast) <- names(fixef(fit))
  A_nams <- paste0(trt_name, A_levels)
  B_nams <- paste0(trt_name, B_levels)
  stopifnot(all(A_nams %in% names(fixef(fit))))
  stopifnot(all(B_nams %in% names(fixef(fit))))
  A_ind <- grep(paste0("(", paste0(A_nams, collapse="|") , ")"), names(fixef(fit)))
  B_ind <- grep(paste0("(", paste0(B_nams, collapse="|") , ")"), names(fixef(fit)))
  diff_contrast[A_ind] <- 1/length(A_levels)
  diff_contrast[B_ind] <- -1/length(B_levels)
  stopifnot(abs(sum(diff_contrast)) < 1e-7)
  diff_contrast
}

get_contrasts_oat <- function(fit){
  contest(fit, rbind(
    `Bison: pure - mixed` = contrast_vector_A_minus_B(fit, A = 4, B=10:12),
    `Lion:  pure - mixed` = contrast_vector_A_minus_B(fit, A = 5, B=13:15),
    `Troll: pure - mixed` = contrast_vector_A_minus_B(fit, A = 6, B=16:18)
  ), joint=FALSE, confint=FALSE)
}

get_contrasts_lupin <- function(fit){
  contest(fit, rbind(
    `Lunabor:  pure - mixed` = contrast_vector_A_minus_B(fit, A = 7, B=c(10,13,16)),
    `Probor :  pure - mixed` = contrast_vector_A_minus_B(fit, A = 8, B=c(11,14,17)),
    `Jowisz :  pure - mixed` = contrast_vector_A_minus_B(fit, A = 9, B=c(12,15,18))
  ), joint=FALSE, confint=FALSE)
}


MODELS <- list(
  M1.o = .aux(
    RHS_formula_core = "is.mixed + cv_oat",
    # RHS_formula_core = "0 + trt",
    analysis_fun = anova
  ),
  M1.l = .aux(
    RHS_formula_core = "is.mixed + cv_lupin",
    # RHS_formula_core = "0 + trt",
    analysis_fun = anova
  ),
  M2.o = .aux(
    # RHS_formula_core = "is.mixed + cv_oat"
    RHS_formula_core = "0 + trt",
    analysis_fun = get_contrasts_oat
  ),
  M2.l = .aux(
    # RHS_formula_core = "is.mixed + cv_lupin"
    RHS_formula_core = "0 + trt",
    analysis_fun = get_contrasts_lupin
  ),
  M3 = .aux(
    RHS_formula_core = "cv_oat * cv_lupin"
  ),
  M4.max.1 = .aux(
    RHS_formula_core = "0 + cv_",
    analysis_fun = pval_FUN_which_coef_diff_from_val(val = 1)
  ),
  M4.max.2 = .aux(
    RHS_formula_core = "0 + cv_",
    analysis_fun = pval_FUN_which_coef_diff_from_val(val = 2)
  ),
  M4.max.01 = .aux(
    RHS_formula_core = "0 + cv_",
    analysis_fun = pval_FUN_which_coef_diff_from_val(val = 0.1)
  ),
  M4.max.09 = .aux(
    RHS_formula_core = "0 + cv_",
    analysis_fun = pval_FUN_which_coef_diff_from_val(val = 0.9)
  ),
  M4.min.0 = .aux(
    RHS_formula_core = "0 + cv_",
    analysis_fun = pval_FUN_which_coef_diff_from_val(val = 0)
  )
)

eval_model <- function(hyp, D){
  if(!is.na(hyp$model_type) & hyp$model_type != "")
    MODELS[[hyp$model_type]](
      hyp, 
      D
      )
}

# eval_model_cache <- function(hyp, D){
#   model_fun <- MODELS[[hyp$model_type]]
#   xfun::cache_rds({
#       eval_model(hyp, D)
#     },
#     hash = list(
#       hyp,
#       D,
#       as.character(functionBody(model_fun)),
#       formals(model_fun),
#       as.character(functionBody(eval_model)),
#       formals(eval_model)
#     ),
#     dir = "cache/",
#     file = paste0("eval_model_", 
#                   hyp$H_ID, "_"
#                   )
#   )
# }

refit_as_lmer <- function(rfit){
  stopifnot(inherits(rfit, "rlmerMod")|inherits(rfit, "lmerMod"))
  suppressMessages(lmer(attr(rfit@frame, "formula"), data=rfit@frame))
}


anova.rlmerMod <- function(rfit, nsim = 20){
  # get null model
  # simulate under null
  # get contrast matrixes
  # get l2 norm of contrasts for beta_star
  # get l2 norm of contrasts for beta_simulated_H0
  # compute p-value

  stop("I am important")
  # refit as lmer
  lfit <- refit_as_lmer(rfit)
  # null with rfit VarCorr
  VarCorr(lfit) <- VarCorr(rfit)
  fixef(lfit)[] <- 0

  # simulate under NULL
  Y_sim <- simulate(lfit, nsim)
  
  # update(lfit, Y_sim[[1]] ~ .)
  Null_fits <- suppressMessages(
    lapply(Y_sim, function(Y_new) update(lfit, Y_new ~ ., data=cbind(rfit@frame, Y_new=Y_new)))
  )
  Null_fixefs <- sapply(Null_fits, fixef)

  # get contrast matrixes
  type3cont <- lmerTest:::get_contrasts_type3(lfit)

  # get l2 norm of contrasts for beta_star
  # get l2 norm of contrasts for beta_simulated_H0
  l2_norms <- matrix(NA, nrow = length(type3cont), ncol = nsim, dimnames = list(names(type3cont), 1:nsim))
  l2_norms_star <- matrix(NA, nrow = length(type3cont), ncol = 5, dimnames = list(names(type3cont), c("star", "H0_mean", "H0_sd", "H0_pval", "H0_pval_fraction")))
  for(cont_name in names(type3cont)){
    l2_norms_star[cont_name, "star"] <- colSums((type3cont[[cont_name]] %*% fixef(rfit))^2)
    l2_norms[cont_name,] <-             colSums((type3cont[[cont_name]] %*% Null_fixefs)^2)
  }
  l2_norms_star[,"H0_mean"] <- rowMeans(l2_norms)
  l2_norms_star[,"H0_sd"] <- apply(l2_norms, 1, sd)
  
  # compute p-value:
  
  l2_norms_star[,"H0_pval_fraction"] <- rowMeans(l2_norms > l2_norms_star[,"star"], na.rm = TRUE)
  # chisq_df <- sapply(type3cont, nrow)
  # beta_sigma <- ???
  # pval that l2_norms_star[,"star"] > l2_norms
  # l2_norms_star[,"H0_pval"] <-  some chi square stuff or figure out distribution of betas: https://math.stackexchange.com/questions/656762/distribution-of-the-sum-of-squared-independent-normal-random-variables

  l2_norms_star |> as.data.frame() |> signif(3)
}


if (FALSE){
  source("R/DataPrep.R")
  source("R/read_hyp.R")
  # attach lmerTest namespace
  library(lmerTest)
  # loadNamespace("lmerTest")
  hyp <- hypotheses["H4.1.l",]
  hyp <- hypotheses["H9",]
  hyp <- hypotheses["H1.l.rob",]
  hyp <- hypotheses["H12",]


  em <- eval_model(hyp, D)
  

  form <- as.formula(paste0(hyp$varname ,"~ 0 + cv_", get_str_formula_rest(D, hyp)))
  lfit <- lmer(form, data=D, subset = is.oat & !is.na(D[[hyp$varname]]), weights=NULL)
  rfit <- fit <- rlmer(form, data=D, subset = is.mixed & !is.na(D[[hyp$varname]]))


  rfit@frame[[1]] <- rfit@frame[[1]] - Dt$diff.b - Dt$diff.e
  update(rfit, rfit@frame[[1]]~., data=rfit@frame)
  update(rfit, weights = robustlmm:::wgt.e(fit))
  lfit2 <- summary(lmer(form, rfit@frame, weights=robustlmm:::wgt.e(fit) * Dt[,"wgt.b"]))
  lfit3 <- summary(update(lfit, competitive_balance_lupins -  - Dt$diff.b - Dt$diff.e ~ .))

  summary(lfit) |> coef()
  summary(rfit) |> coef()
  summary(lfit2) |> coef() 
  summary(lfit3) |> coef() 

  coef(summary(lfit))
  lmerTest:::contest1D.lmerModLmerTest()
  COEFs <- coef(summary(rfit))
  pval_FUN_which_coef_diff_from_val()(coef(summary(rfit)))

  rfit@theta
  as_lmerModLmerTest
  lmerTest:::as_lmerModLT(rfit)
lmerTest:::as_lmerModLT
lmerTest:::get_covbeta

  nullfit <- lfit
  fixef(nullfit)[] <- 0
  Y_sim <- simulate(nullfit, 1, newdata = rfit@frame)
  tmp <- cbind(rfit@frame, nully=Y_sim[[1]])
  nullfit <- update(nullfit, nully ~ ., data=tmp)
  rnullfit <- rlmer(attr(nullfit@frame, "formula"), data=tmp)
  cbind(anova(rnullfit, 1)[,"H0_pval"]> anova(nullfit)[,"Pr(>F)"])
  

  hyp <- hypotheses["H2.o",]
  hyp <- hypotheses["H4.1.l",]
  hyp <- hypotheses["H12",]
  source("R//ggplot_hyp.R")
  
  library(censReg)
  censReg(formula, left = min(y), right = max(y), data = ,
    subset = NULL, logLikOnly = FALSE, ... )

  data("Gasoline", package = "plm")
  library(plm)
    Gas <- pdata.frame(Gasoline, index = c("country", "year"), drop.index = TRUE)

  names(D)
  D$date_when_covered |> is.na()
  D$date_when_covered_diff_mix_to_pure_lupin |> is.na()
  D$date_when_covered_diff_mix_to_pure_oat |> is.na()
  d <- D$date_when_covered
  summary(d)
  D$day_of_year_when_covered <- as.integer(strftime(D$date_when_covered, format = "%j"))
  # plot(a ~ year + season + location, data=D)
  # -> plot seperately for seasons
  
  

  huberM(na.omit(c(1:100, 10^9)), se = TRUE)
  HuberM(na.omit(c(1:100, 10^9)), conf.level = 0.95)


  y <- na.omit(y)
  HuberM(y, conf.level = 0.95, na.rm = TRUE)
  HuberM(na.omit(y))
  Hub <- huberM(na.omit(y), se=TRUE)
  Hub$mu - c(-1, 1)*Hub$SE*qt(0.975, length(na.omit(y)-1))
  HuberM(na.omit(y), conf.level = 0.95, ci.type = "boot")

  get_man_ci_huber <- function(y, conf.level = 0.95, s=mad(y, na.rm=TRUE)){
    y <- na.omit(y)
    Hub <- huberM(y, se=TRUE, s=s)
    c(
      mu = Hub$mu,
      ci = Hub$mu - c(-1, 1)*Hub$SE*qt(conf.level/2, length(y)-1)
    )
  }


  aggregate(date_when_covered_diff_mix_to_pure_lupin ~ trt, data=D, FUN=HuberM, conf.level = 0.95, s = mad(D$date_when_covered_diff_mix_to_pure_lupin, na.rm=TRUE), ci.type="boot")


  aggregate(date_when_covered_diff_mix_to_pure_lupin ~ trt, data=D, FUN=get_man_ci_huber, conf.level = 0.95, s = mad(D$date_when_covered_diff_mix_to_pure_lupin, na.rm=TRUE))

  table(D$date_when_covered_diff_mix_to_pure_lupin)
  median(y,na.rm=T)
  MedianCI(y, na.rm=T)
}


# Custom print method
print.model.analysis <- function(x){
  is.knitr <- !is.null(getOption("knitr.in.progress"))

  if(is.knitr) cat("\n\n```{r}\n")

  cat("VarName: ", x$varname, "\n")

  # Formula
  if (length(x$formula == 3)){
    print(x$formula[c(1,3)], showEnv=FALSE)
  } else {
    print(x$formula, showEnv=FALSE)
  }

  # Results
  cat("\n")
  print(x$results, signif.legend=FALSE)

  # p-values
  if (!is.null(x$p_vals)){
    cat("\n")
    cat("P-values: \n")
    print(x$p_vals)
  }

  # emmeans posthoc
  if (!is.null(x$emmeans_posthoc)) {
    cat("\n\n**Post-hoc pairwise comparisons (emmeans):**\n\n")
    print(x$emmeans_posthoc)
  }

  if(is.knitr) cat("\n```\n")
}



######################################################
# TESTs-----------------------------------------------
if(FALSE){
  source("R/DataPrep.R")
  source("R/read_hyp.R")
  eval_model(hypotheses["H4.1.l",], D)
  hypotheses$model_type
  for (i in 1:nrow(hypotheses)){
    print(i)
    Sys.sleep(1)
    hyp <- hypotheses[i,]
    print(eval_model(hyp, D))
  }


  # TEST CONTRASTS
    # OAT : pure - mixed
    hyp <- hypotheses["H4.1.o",]
    ff <- get_str_formula_rest(D, hyp)
    fit <- lmer(as.formula(paste(hyp$varname, " ~ 0 + trt", ff)), data=D, subset=is.oat)
    contest(fit, rbind(
      `Bison: pure - mixed` = contrast_vector_A_minus_B(fixef(fit), A = 4, B=11:14),
      `Lion:  pure - mixed` = contrast_vector_A_minus_B(fixef(fit), A = 5, B=15:18),
      `Troll: pure - mixed` = contrast_vector_A_minus_B(fixef(fit), A = 5, B=19:22)
    ), joint=FALSE)
    get_contrasts_oat(fit)
    # lupin : pure - mixed
    hyp <- hypotheses["H4.1.l",]
    ff <- get_str_formula_rest(D, hyp)
    fit <- lmer(as.formula(paste(hyp$varname, " ~ 0 + trt", ff)), data=D, subset=is.lupin)
    get_contrasts_lupin(fit)
    contest(fit, rbind(
      `Anicia:  pure - mixed` = contrast_vector_A_minus_B(fixef(fit), A = 7, B=c(11,15,19)),
      `Beluga:  pure - mixed` = contrast_vector_A_minus_B(fixef(fit), A = 8, B=c(12,16,20)),
      `Spaeths: pure - mixed` = contrast_vector_A_minus_B(fixef(fit), A = 9, B=c(13,17,21)),
      `Itaca:   pure - mixed` = contrast_vector_A_minus_B(fixef(fit), A = 10, B=c(14,18,22))
    ), joint=FALSE)

  # sanity check coefficient interpretations  
    hyp <- hypotheses["H4.1.l",]
    ff <- get_str_formula_rest(D, hyp)
    fit <- lmer(as.formula(paste(hyp$varname, " ~ 0 + trt", ff)), data=D, subset=is.oat)
    coef(summary(fit))
    aggregate(D[[hyp$varname]] ~ trt, data=D, FUN=mean, subset=is.oat)
    aggregate(D[[hyp$varname]] ~ year, data=D, FUN=mean, subset=is.oat)
    aggregate(D[[hyp$varname]] ~ location, data=D, FUN=mean, subset=is.oat)
}
