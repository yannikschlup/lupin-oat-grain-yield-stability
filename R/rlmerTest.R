# Author: Lukas Graz             Date: 2024-11-14
# -----------------------------------------------
# P-Values for rlmer via lmerTest
# -----------------------------------------------

#' How it works
#' ============
#' 
#' lmerTest computes ALL p-values via conrast1D and contestMD functions.
#' those functions in return only use:
#' - beta
#' - Jac_list     = refit_as_lmer(rlmer_object)@Jac_list
#' - vcov_varpar  = refit_as_lmer(rlmer_object)@vcov_varpar
#' - vcov_beta    = vcov(rlmer_object)
#' - sigma        = sigma(rlmer_object)
#' We define a new class, rlmerModTest, that contains those slots. and are 
#' obtained as indicated above.

library(robustlmm)
library(lmerTest)
library(simr)

setClass(
  "rlmerModTest",
  contains="rlmerMod",
  slots=c(
    Jac_list = "list",
    vcov_varpar = "matrix",
    vcov_beta = "matrix",
    sigma = "numeric"
    )
) -> rlmerModTest


refit_as_lmer <- function(object){
  if (inherits(object, "lmerMod"))
    return(object)
  stopifnot(inherits(object, "rlmerMod"))
  suppressMessages(lmer_refit <- lmerTest::lmer(attr(object@frame, "formula"), data=object@frame))
  lmer_refit@vcov_varpar
  lmer_refit@theta
  object@theta
  lmer_refit@vcov_beta <- as.matrix(vcov(object))
}



contestMD.rlmerModTest <- lmerTest:::contestMD.lmerModLmerTest


contest1D.rlmerModTest <- lmerTest:::contest1D.lmerModLmerTest


# anova.rlmerMod <- function (object, type = c("III", "II", "I", "3", "2", "1"), 
#     ddf = c("Satterthwaite", "Kenward-Roger", "lme4"), silent=TRUE) {
#     stopifnot(inherits(object, "rlmerMod"))
#     ddf <- match.arg(ddf)
#     object <- as.rlmerModTest(object)
#     if (silent){
#       suppressWarnings(
#         suppressMessages(
#           return(
#             lmerTest:::single_anova(object = object, type = type, ddf = ddf)
#     )))} 
#     lmerTest:::single_anova(object = object, type = type, ddf = ddf)
# }


summary.rlmerMod <- function(object){
    resp <- object@resp
    devC <- object@devcomp
    dd <- devC$dims
    cmp <- devC$cmp
    sig <- sigma(object)
    coefs <- cbind(Estimate = fixef(object), `Std. Error` = sig * 
        sqrt(diag(object@pp$unsc())))
    if (nrow(coefs) > 0) {
        coefs <- cbind(coefs, df = lmerTest:::get_coefmat(refit_as_lmer(object))[,"df"]) # new
        coefs <- cbind(coefs, `t value`=coefs[, 1]/coefs[, 2], deparse.level = 0)
        coefs <- cbind(coefs, `Pr(>|t|)` = 2 * pt(-abs(coefs[, 
            "t value"]), df = coefs[, "df"]))
    }
    varcor <- VarCorr(object)
    structure(list(methTitle = robustlmm:::.methTitle(object), devcomp = devC, 
        ngrps = sapply(object@flist, function(x) length(levels(x))), 
        coefficients = coefs, sigma = sig, vcov = vcov(object, 
            correlation = TRUE, sigm = sig), varcor = varcor, 
        call = object@call, wgt.e = robustlmm:::wgt.e(object), wgt.b = robustlmm:::wgt.b(object), 
        rho.e = robustlmm:::rho.e(object), rho.sigma.e = robustlmm:::rho.e(object, "sigma"), 
        rho.b = robustlmm:::rho.b(object), rho.sigma.b = robustlmm:::rho.b(object, "sigma"), 
        residuals = residuals(object, scaled = TRUE)), class = "summary.rlmerMod")
}

# # remove inheritance check from lmerTest:::get_model_matrix
# get_model_matrix <- lmerTest:::get_model_matrix
# body(get_model_matrix)[[3]] <- NULL
# # assign get_model_matrix to lmerTest namespace
# unlockBinding("get_model_matrix", getNamespace("lmerTest"))
# assignInNamespace("get_model_matrix", get_model_matrix, ns = "lmerTest")
# lockBinding("get_model_matrix", getNamespace("lmerTest"))
# rm(get_model_matrix)



as.rlmerModTest <- function(object, ...){
  if(inherits(object, "rlmerModTest"))
    return(object)
  if(inherits(object, "rlmerMod")){
    obj <- as(object, "rlmerModTest")
    refit_lmerModLmerTest <- refit_as_lmer(object)
    obj@Jac_list <- refit_lmerModLmerTest@Jac_list
    obj@vcov_varpar <- refit_lmerModLmerTest@vcov_varpar
    obj@vcov_beta <- as.matrix(vcov(object))
    obj@sigma <- sigma(object)
  } else {
    stop("object must be of class rlmerMod or rlmerModTest")
  }
  # as_dummy(obj)
  obj
}




# as_dummy <- function (model, tol = 1e-08) {
#     # if (!inherits(model, "lmerMod")) 
#     #     stop("model not of class 'lmerMod': cannot coerce to class 'lmerModLmerTest")
#     mc <- getCall(model)
#     args <- c(as.list(mc), devFunOnly = TRUE)
#     if (!"control" %in% names(as.list(mc))) 
#         args$control <- lme4::lmerControl(check.rankX = "silent.drop.cols")
#     Call <- as.call(c(list(quote(lme4::lmer)), args[-1]))
#     ff <- environment(formula(model))
#     pf <- parent.frame()
#     sf <- sys.frames()[[1]]
#     ff2 <- environment(model)
#     devfun <- tryCatch(eval(Call, envir = pf), error = function(e) {
#         tryCatch(eval(Call, envir = ff), error = function(e) {
#             tryCatch(eval(Call, envir = ff2), error = function(e) {
#                 tryCatch(eval(Call, envir = sf), error = function(e) {
#                   "error"
#                 })
#             })
#         })
#     })
#     if ((is.character(devfun) && devfun == "error") || !is.function(devfun) || 
#         names(formals(devfun)[1]) != "theta") 
#         stop("Unable to extract deviance function from model fit")
#     dummy(model, devfun, tol = tol)
# }



# dummy <- function (model, devfun=NULL, tol = 1e-08) {
#     is_reml <- getME(model, "is_REML")
#     # res <- as(model, "lmerModLmerTest")
#     res <- as(model, "rlmerModTest") ### new
#     res@sigma <- sigma(model)
#     res@vcov_beta <- as.matrix(vcov(model))
#     varpar_opt <- unname(c(res@theta, res@sigma))
#     h <- numDeriv:::hessian(func = pbkrtest:::devfun_vp, x = varpar_opt, 
#         devfun = devfun, reml = is_reml)
#     eig_h <- eigen(h, symmetric = TRUE)
#     evals <- eig_h$values
#     neg <- evals < -tol
#     pos <- evals > tol
#     zero <- evals > -tol & evals < tol
#     if (sum(neg) > 0) {
#         eval_chr <- if (sum(neg) > 1) 
#             "eigenvalues"
#         else "eigenvalue"
#         evals_num <- paste(sprintf("%1.1e", evals[neg]), collapse = " ")
#         warning(sprintf("Model failed to converge with %d negative %s: %s", 
#             sum(neg), eval_chr, evals_num), call. = FALSE)
#     }
#     if (sum(zero) > 0) {
#         eval_chr <- if (sum(zero) > 1) 
#             "eigenvalues"
#         else "eigenvalue"
#         evals_num <- paste(sprintf("%1.1e", evals[zero]), collapse = " ")
#         warning(sprintf("Model may not have converged with %d %s close to zero: %s", 
#             sum(zero), eval_chr, evals_num))
#     }
#     pos <- eig_h$values > tol
#     q <- sum(pos)
#     h_inv <- with(eig_h, {
#         vectors[, pos, drop = FALSE] %*% diag(1/values[pos], 
#             nrow = q) %*% t(vectors[, pos, drop = FALSE])
#     })
#     res@vcov_varpar <- 2 * h_inv
#     ## NO LOGLIK:
#     # Jac <- numDeriv::jacobian(func = lmerTest:::get_covbeta, x = varpar_opt, 
#     #     devfun = devfun)
#     # res@Jac_list <- lapply(1:ncol(Jac), function(i) array(Jac[, 
#     #     i], dim = rep(length(res@beta), 2)))
#     res
# }



###################################################
# TESTS
###################################################
if (FALSE){
  data(sleepstudy)
  sleepstudy$FDays <- factor(sleepstudy$Days)

  lfit <- lmer(Reaction ~ FDays + (1 | Subject), data=sleepstudy)
  rfit <-rlmer(Reaction ~ FDays + (1 | Subject), data=sleepstudy)

  coef(summary(lfit))
  coef(summary(rfit))

  get_last_column_from_df(anova(rfit))
  get_last_column_from_df(anova(lfit))

  

  # simulate under null-data (using sleepstudy as basis) 
  library(simr)
  nsim <- 50
  nullfit <- lfit
  fixef(nullfit)[] <- 0
  set.seed(123)
  Y_sim <- simulate(nullfit, nsim)

  # fit lmer and rlmer
  l_sim <- lapply(Y_sim, function(Y_new) 
    lmer(Y_new ~ FDays + (1 | Subject), data=cbind(sleepstudy, Y_new=Y_new)))

  r_sim <- lapply(Y_sim, function(Y_new) 
    rlmer(Y_new ~ FDays + (1 | Subject), data=cbind(sleepstudy, Y_new=Y_new)))

  # sapply(l_sim, fixef)|> t() |> summary()
  # sapply(r_sim, fixef)|> t() |> summary()

  # get p-values
  r_pvals <- sapply(r_sim, \(x) anova(x)[,"Pr(>F)",drop=TRUE])
  l_pvals <- sapply(l_sim, \(x) anova(x)[,"Pr(>F)",drop=TRUE]) 

  ks.test(r_pvals,"punif",0,1)
  ks.test(l_pvals,"punif",0,1)


  hist(r_pvals)
  hist(l_pvals)

  # compare p-values
  cor(l_pvals, r_pvals)



  # ---------------------
  # Do the same again, but take a non-balanced subset of the data
  sleepstudy2 <- sleepstudy[sample(nrow(sleepstudy), nrow(sleepstudy), replace=TRUE),]
  lfit <- lmer(Reaction ~ FDays + (1 | Subject), data=sleepstudy2)
  nullfit <- lfit
  fixef(nullfit)[] <- 0
  set.seed(123)
  Y_sim <- simulate(nullfit, nsim)

  # fit lmer and rlmer
  l_sim <- lapply(Y_sim, function(Y_new) 
    lmer(Y_new ~ FDays + (1 | Subject), data=cbind(sleepstudy2, Y_new=Y_new)))

  r_sim <- lapply(Y_sim, function(Y_new) 
    rlmer(Y_new ~ FDays + (1 | Subject), data=cbind(sleepstudy2, Y_new=Y_new)))

  # sapply(l_sim, fixef)|> t() |> summary()
  # sapply(r_sim, fixef)|> t() |> summary()

  # get p-values
  r_pvals <- sapply(r_sim, \(x) anova(x)[,"Pr(>F)",drop=TRUE])
  l_pvals <- sapply(l_sim, \(x) anova(x)[,"Pr(>F)",drop=TRUE]) 

  ks.test(r_pvals,"punif",0,1)
  ks.test(l_pvals,"punif",0,1)

  # compare p-values
  cor(l_pvals, r_pvals)
}
