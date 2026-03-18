Sys.setenv(RENV_PATHS_PROJECT = "../")
renv::activate()
#source("renv-path.R")

# # Parallelization
# library(parallel)
# detectCores()
# cl <- makeCluster(4)
# #run calculations
# stopCluster(cl)
