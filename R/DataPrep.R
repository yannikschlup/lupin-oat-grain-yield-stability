latest_file <- "data/251021_MastersheetLupinOatNnoN.xlsx"

library(readxl)
# DATA
D_orig <- read_xlsx(
  latest_file, 
  sheet = 1, 
  skip = 11, 
  col_names = TRUE,
  na = c("", "NA", "no", "No", "#VALUE!")
  )

# METADATA
D_meta <- read_xlsx(latest_file, sheet = 1, n_max = 11, col_names = names(D_orig)) |> as.matrix()
rownames(D_meta) <- D_meta[,"year"]
D_meta[,"year"] <- NA
D_meta["Data type", "year"] <- "factor"
D_meta <- as.data.frame(t(D_meta))
D_meta$`Data type` |> table()
D_meta$`Data type`[D_meta$`Data type` %in% c(
  "bbch difference",
  "bbch score",
  "days difference",
  "numeric",
  "percent",
  "percentage"
)] <- "numeric"


### TYPE ENCODING
D_meta[, "Data type"] |> table() |> cbind()
D_meta$unique_values <- sapply(D_orig, \(x) length(unique(x)))
# View(D_meta[,c("Data type", "unique_values")])

# define types
factors_cols <- names(D_orig)[D_meta$`Data type` %in% c("factor")]
num_cols <- names(D_orig)[D_meta$`Data type` %in% c("numeric", "fraction", "days difference", "percent", "ordinal", "likert")]
date_cols <- names(D_orig)[D_meta$`Data type` %in% c("date")]
stopifnot(length(factors_cols) + length(num_cols) + length(date_cols) == ncol(D_orig))

D_typed <- D_orig[!(D_orig$treatment %in% 19:20), ]

# encode types correctly:
D_typed[factors_cols] <- lapply(D_typed[factors_cols], as.factor)
D_typed[num_cols] <- lapply(D_typed[num_cols], as.numeric)
# D_typed[date_cols] <- lapply(D_typed[date_cols], as.Date)

### TYPE ENCODING END

# # set +/- 1e10 to +/- 50
# for (var in c(
#   "date_when_covered_diff_mix_to_pure_oat", 
#   "date_when_covered_diff_mix_to_pure_lupin"
# )){
#   ind <- !is.na(D_typed[[var]]) & abs(D_typed[[var]]) == 1e10
#   D_typed[[var]][ind] <- sign(D_typed[[var]][ind]) * 50
# }

# # add intercept for weed_vol_percent   0.5%
# stopifnot(all(D_typed$weed_vol_percent_june >= 0, na.rm = TRUE))
# D_typed$weed_vol_percent_june <- D_typed$weed_vol_percent_june + 0.00005
# D_typed$lupin_athracnose_no_pods <- D_typed$lupin_athracnose_no_pods + 0.2
# # add day of year
# D_typed$day_of_year_when_covered <- as.integer(strftime(D_typed$date_when_covered, format = "%j"))
# # set NA's to 200
# D_typed$day_of_year_when_covered[is.na(D_typed$day_of_year_when_covered)] <- 200

D_typed$cv_ <- interaction(D_typed$cv_oat, D_typed$cv_lupin, drop = TRUE)

D_typed$trt <- D_typed$treatment

# alias season, year, location with s y l
D_typed$y <- D_typed$year
D_typed$l <- D_typed$location

D <- D_typed
D

### is.???
is.pure_oat_nF<-D$is.pure_oat_nF <-D$treatment %in% 1:3  # no fertilizer
is.pure_oat   <-D$is.pure_oat    <-D$treatment %in% 4:6 # fertiliser
is.pure_lupin <-D$is.pure_lupin  <-D$treatment %in% 7:9
is.mixed      <-D$is.mixed       <-D$treatment %in% 10:18
is.pure       <-D$is.pure        <-is.pure_oat | is.pure_lupin
is.oat        <-D$is.oat         <-is.pure_oat | is.mixed # no unfertilized pure oat
is.lupin      <-D$is.lupin       <-is.pure_lupin | is.mixed
is.all        <-D$is.all         <-is.pure_oat | is.pure_lupin | is.mixed
