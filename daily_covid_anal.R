## download daily update COVID-19 data from data.go.kr
## and add it to the existing data files (dat.csv, dat.rds)
# eval(parse('covid_data_daily_update.R', 'UTF-8'))
library(XML)
library(RCurl)
library(dplyr)
# Sys.setlocale("LC_ALL", "Korean")
# options(encoding = "UTF-8")
## data portal
service_url <- "http://openapi.data.go.kr/openapi/service/rest/Covid19/getCovid19SidoInfStateJson"
rows <- 4 #number of rows
## page
pg <- 1
# service key (register at https://www.data.go.kr/data/15043376/openapi.do)
# web site also provides a api request tool to check if the service key works
# decoded key works for the web site but in R encoded key (with the percent sign works)
# service_key <- "KNu/yp3/Km2MWttRvbZXNE+b/EH44cfIYDvMgPFoSKUuKZIxqQoGu0gnAtwslGjupx+E6vp0bwf/SRPcLfXSYQ=="
service_key <- "KNu%2Fyp3%2FKm2MWttRvbZXNE%2Bb%2FEH44cfIYDvMgPFoSKUuKZIxqQoGu0gnAtwslGjupx%2BE6vp0bwf%2FSRPcLfXSYQ%3D%3D"
## service period
start_dt <- gsub("-", "", as.character(Sys.Date()))
end_dt <- start_dt

# In case the program didn't run for any reason, you can set the start_dt and
# end_dt manually, and run the program.
# Also, you have to ensure that the last row of the data_full has the date one
# day before the start_dt

# library(tidyverse)
# readRDS("outputs/data_full.rds") %>%
#   dplyr::filter(createDt < as.Date("2021-07-07")) %>%
#   saveRDS("outputs/data_full.rds")
#
# start_dt <- "20220116"
# end_dt <- "20220117"

uri <-  paste0(service_url,
               paste0("?serviceKey=", service_key),
               paste0("&pageNo=", pg),
               paste0("&numOfRows=", rows),
               paste0("&startCreateDt=", start_dt),
               paste0("&endCreateDt=", end_dt))

xml_doc <- xmlTreeParse(uri, useInternalNodes = TRUE, encoding = "UTF-8")
root_node <- xmlRoot(xml_doc)
xml_data <- xmlToDataFrame(nodes = getNodeSet(root_node, '//item'))

# xml_data <- xml_data[order(xml_data$stdDay), ]
# saveRDS(ds, "outputs/data_full_20200302_20200531.rds")

# for(i in 1:(nrow(dd_total)-1)) {
#   if(as.numeric(as.Date(dd_total[i+1, 1]) - as.Date(dd_total[i, 1])) != 1) {
#     message(paste0("i = ", i))
#   }
# }
# which(as.Date(dd$createDt) == "2020-10-12")
# which(as.Date(dd$createDt) == "2020-10-13")
#
# rr <- dd[2564, ]
# for(i in c(2:3, 7:11)) dd[, i] <- as.numeric(dd[, i])
# rr[, c(2:3, 7:11)] <- colSums(dd[2528:2545, c(2:3, 7:11)])
# rr[, 12:13] <- NA
# rr[, c(1,14)] <- dd[2545, c(1,14)]
# dd <- rbind(dd[1:2545, ], rr, dd[2546:nrow(dd), ])
# saveRDS(dd, "outputs/data_full_20200601_20201231.rds")

# data_full <- readRDS("outputs/data_full.rds")
# data_full <- data_full[data_full$createDt < as.Date("2022-01-16"),]

ncat <- 19 # 19 categories for regions: 17 regions, imported, and total
days <- as.numeric(as.Date(end_dt, "%Y%m%d") - as.Date(start_dt, "%Y%m%d") + 1)
if (nrow(xml_data) == ncat * days) {
  ## original data with full variables
  data_full <- readRDS("outputs/data_full.rds")
  # saveRDS(data_full, "outputs/data_full_3Dec2021.rds")
  xml_data <- xml_data[order(xml_data$stdDay), ] # make ascending date
  if (as.Date(tail(xml_data$createDt, 1)) > as.Date(tail(data_full$createDt, 1))) {
    data_full <- rbind(data_full, xml_data)
    saveRDS(data_full, "outputs/data_full.rds")
  }
}


# data_full <- readRDS("outputs/data_full.rds")
dat <- data_full %>% filter(gubunEn == "Total") %>%
  select(createDt, localOccCnt)
names(dat) <- c("date", "daily_confirmed")
dat$date <- as.Date(dat$date)
dat$daily_confirmed <- as.numeric(dat$daily_confirmed) # originally character
saveRDS(dat, "outputs/dat.rds")
## save the file also in csv format for easier handling in the shiny app
readr::write_csv(dat, "daily_sim/dat.csv")

# ## Gyeonggi-do
# dat_gg <- data_full %>% filter(gubunEn == "Gyeonggi-do") %>%
#   select(createDt, localOccCnt)
# names(dat_gg) <- c("date", "daily_confirmed")
# dat_gg$date <- as.Date(dat_gg$date)
# dat_gg$daily_confirmed <- as.numeric(dat_gg$daily_confirmed)
# saveRDS(dat_gg, "daily_sim/dat_gg.rds")
# ## save the file also in csv format for easier handling in the shiny app
# data.table::fwrite(dat_gg, "daily_sim/dat_gg.csv")

# devtools::install()
library(pfilter)
library(parallel)
library(doParallel)
library(foreach)

# latent variables for 1 Jan 2021 estimated by  previous particle filtering
# simulations that used data from 1 Jun 2021
y0_20210531 <- readRDS("outputs/y0_20210531.rds")
y0 <- round(y0_20210531)
# dat <- readRDS("daily_sim/dat.rds")
# theta[["betavol"]] <- 0.1
# theta[["gamma"]] <- 1/2.5
# usethis::use_data(theta, overwrite = T)
# devtools::load_all(".")

set.seed(23)

ncores <- detectCores()
cl <- makeCluster(getOption("cl.cores", ncores - 2))
doParallel::registerDoParallel(cl)

npart <- 1e4
nrep <- 1e3
dt <- 0.2

pf <- foreach(i = 1:nrep, .packages = "pfilter", .inorder = F) %dopar% {
  extract_trace(params = theta,
                y = y0,
                data = dat,
                data_type = "confirmation",
                npart = npart,
                tend = nrow(dat),
                dt = dt,
                error_pdf = "negbin",
                negbin_size = 30,
                stoch = FALSE)
}
parallel::stopCluster(cl)

# saveRDS(pf, "daily_sim/pf.rds")
saveRDS(pf, "outputs/pf.rds")
Rt <- as.data.frame(sapply(pf, function(x) x[, "Rt"]))
daily_conf <- as.data.frame(sapply(pf, function(x) x[, "CR"]))
## save as the csv files as they are easier to
# data.table::fwrite(Rt, "daily_sim/Rt.csv")
# data.table::fwrite(daily_conf, "daily_sim/daily_confirmed.csv")

set.seed(23)
nsample <- 200
ids <- sample.int(nrep, nsample)

nstates <- 9 # SEPIR, CE, CI, CR, A
str <- names(pf[[1]])
samp <- lapply(str[1:nstates], function(x) extract_sample_pf(pf, x, days = 1, ids))
smpl_last_states <- as.data.frame(do.call(cbind, samp))
names(smpl_last_states) <- str[1:nstates]

# samp_Rt <- sapply(str[nstates + 1], function(x) extract_sample_pf(pf, x, days, ids))
# smpl_Rt <- as.data.frame(samp_Rt)
smpl_Rt <- Rt[, ids]
smpl_daily_confirmed <- daily_conf[, ids]

# sub_pf <- pf[ids]
# smpl_daily_confirmed <- as.data.frame(sapply(sub_pf, function(x) x[, "CR"]))

data.table::fwrite(smpl_last_states, "daily_sim/smpl_last_states.csv")
data.table::fwrite(smpl_Rt, "daily_sim/smpl_Rt.csv")
data.table::fwrite(smpl_daily_confirmed, "daily_sim/smpl_daily_confirmed.csv")

## Gyeonggi-do
# dat_gg <- readRDS("daily_sim/dat_gg.rds")
# y20210101_gg <- readRDS("outputs/y20210101_gg.rds")
# ncores <- detectCores()
# cl <- makeCluster(getOption("cl.cores", ncores - 2))
# doParallel::registerDoParallel(cl)
#
# pf_gg <- foreach(i = 1:1e3, .packages = "pfilter", .inorder = F) %dopar% {
#   extract_trace(params = theta,
#                 y = y20210101_gg,
#                 data = dat_gg,
#                 data_type = "confirmation",
#                 npart = 1e4,
#                 tend = nrow(dat_gg),
#                 dt = 0.25,
#                 error_pdf = "negbin",
#                 negbin_size = 30,
#                 stoch = TRUE)
# }
# parallel::stopCluster(cl)
# saveRDS(pf_gg, "daily_sim/pf_gg.rds")
# Rt_gg <- as.data.frame(sapply(pf_gg, function(x) x[, "Rt"]))
# daily_conf_gg <- as.data.frame(sapply(pf_gg, function(x) x[, "CR"]))
# ## save as the csv files as they are easier to
# data.table::fwrite(Rt_gg, "daily_sim/Rt_gg.csv")
# data.table::fwrite(daily_conf_gg, "daily_sim/daily_confirmed_gg.csv")

