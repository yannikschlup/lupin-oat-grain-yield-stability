# Author: Lukas Graz             Date: 2024-10-03
# -----------------------------------------------
# Check if there were any errors in excel handeling (which would introduce recurring sequences like those detected aboves)
# -----------------------------------------------


source("./R/DataPrep.R")

# GET RECURRING SEQUENCES:
recurring <- xfun::cache_rds({
    get_recurring_seq <- function(v, k=10){
      DF <- NULL
      count <- 0
      for(i in 1:(length(v)-k)){
        for(j in (i+1):(length(v)-k+1)){
          first <- v[i:(i+k-1)]
          second <- v[j:(j+k-1)]
          if(identical(first,second)){
            if (length(unique(first)) == 1){ # skip seq like: 1 1 1 1
              next
            }
            DF <- rbind(DF, data.frame(
              first = i,
              second = j,
              seq1 = paste(first , collapse = ","),
              seq2 = paste(second, collapse = ",")
              )
            )
            count <- count + 1
            # if(count > 20) return(DF)
          }
        }
      }
      return(DF)
    }
    lapply(D_orig, get_recurring_seq)
  },
  file="recurring_seq_list.rds",
  dir = "cache/",
  hash = list(D_orig)
  )


# hier nur die Sequenzen, die in mit dem gleichen treatment starten
recurring2 <- lapply(recurring, function(df) {
  if(is.null(df)) return(NULL)
  xx<<- df
  first_trt  <- D_orig[df$first , "treatment"][[1]]
  second_trt <- D_orig[df$second, "treatment"][[1]]
  df$first_trt <- first_trt
  df$second_trt <- second_trt
  ind <- (first_trt == second_trt) &
         (first_trt %in% c(1,11)) 
  df[ind, ]
})
lapply(recurring2, head)
sapply(recurring2, nrow) |> unlist() |> cbind()


### RECURRING END
