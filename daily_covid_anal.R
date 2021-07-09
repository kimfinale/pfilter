# # download the data from the coronaboard, which is deaily updated
# d <- readr::read_csv("https://raw.githubusercontent.com/jooeungen/coronaboard_kr/master/kr_regional_daily.csv")
# ## get the number of locally transmitted cases (excluding "검역")
# library(dplyr)
# d %>% filter(region != "검역") %>%
#   group_by(as.factor(date)) %>%
#   summarize(confirmed = sum(confirmed)) -> d2
# names(d2) <- c("date", "cumul_confirmed")
# d2$date <- as.Date(as.character(d2$date), format = "%Y%m%d")
# d2 %>% filter(date >= as.Date("2020-12-31")) -> d3
# dat <- d3 %>% filter(date >= as.Date("2021-01-01"))
# dat$daily_confirmed <- diff(d3$cumul_confirmed)
# saveRDS(dat, "daily_sim/dat.rds")
# ## save the file also in csv format for easier handling in the shiny app
# readr::write_csv(dat, "daily_sim/dat.csv")

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

# start_dt <- "20210626"
# end_dt <- "20210626"

uri <-  paste0(service_url,
               paste0("?serviceKey=", service_key),
               paste0("&pageNo=", pg),
               paste0("&numOfRows=", rows),
               paste0("&startCreateDt=", start_dt),
               paste0("&endCreateDt=", end_dt))

xml_doc <- xmlTreeParse(uri, useInternalNodes = TRUE, encoding = "UTF-8")
root_node <- xmlRoot(xml_doc)
xml_data <- xmlToDataFrame(nodes = getNodeSet(root_node, '//item'))

ncat <- 19 # 19 categories for regions: 17 regions, imported, and total
days <- as.numeric(end_dt) - as.numeric(start_dt) + 1
if (nrow(xml_data) == (ncat * days)) {
  ## original data with full variables
  data_full <- readRDS("outputs/data_full.rds")
  xml_data <- xml_data[order(xml_data$stdDay), ] # make ascending date
  if (as.Date(tail(xml_data$createDt, 1)) > as.Date(tail(data_full$createDt, 1))) {
    data_full <- rbind(data_full, xml_data)
    saveRDS(data_full, "outputs/data_full.rds")
  }
}


# check the data that is just created if they are in the right format
# data_full %>% group_by(as.Date(createDt)) %>% summarise(row = n())


  ## data for particle filter fitting (date and daily confirmed)
  # d <- dplyr::select(xml_data, createDt, gubunEn, localOccCnt)
  # d$createDt <- as.Date(d$createDt)
  # d$localOccCnt <- as.numeric(d$localOccCnt)
  # if (tail(d$gubunEn, 1) == "Total") {
  #   # created date, cumulative confirmed cases, locally transmitted cases on the day
  #   d <- d[nrow(d), c("createDt", "localOccCnt")]
  # }

# temporary dat creation from data_full (18Apr2021)
# realized that cases from each region maybe include imported cases
# Thsi is not obvious from daily reports from KDCA, but can be extracted using data.go.kr api
# dat <- data_full %>% filter(gubunEn == "Total") %>%
#   select(createDt, localOccCnt)
# names(dat) <- c("date", "daily_confirmed")
# dat$date <- as.Date(dat$date)
# saveRDS(dat, "daily_sim/dat.rds")

dat <- data_full %>% filter(gubunEn == "Total") %>%
  select(createDt, localOccCnt)
names(dat) <- c("date", "daily_confirmed")
dat$date <- as.Date(dat$date)
dat$daily_confirmed <- as.numeric(dat$daily_confirmed) # orignally character

saveRDS(dat, "daily_sim/dat.rds")
## save the file also in csv format for easier handling in the shiny app
readr::write_csv(dat, "daily_sim/dat.csv")

# dat <- readRDS("daily_sim/dat.rds") # last
# names(d) <- names(dat)

# if (exists(d)) {
#   ## daily_confirmed = locally transmitted cases
#   ## cumul_confirmed = local + imported cases
#   dat_added <- FALSE
#   if (tail(dat$date, 1) < d$date) {
#     dat <- rbind(dat, d)
#     dat_added <- TRUE
#   }
#   if (dat_added) {
#     saveRDS(dat, "daily_sim/dat.rds")
#     ## save the file also in csv format for easier handling in the shiny app
#     readr::write_csv(dat, "daily_sim/dat.csv")
#   }
#

# ## Fit starting from January 1, 2022
library(pfilter)
devtools::load_all(".")
library(parallel)
library(doParallel)
library(foreach)
# latent variables for 1 Jan 2021 estimated by  previous particle filtering
# simulations that used data from 20 Jan 2020
y20210101 <- readRDS("outputs/y20210101.rds")
theta["betavol"] <- 0.10
usethis::use_data(theta, overwrite = T)
devtools::load_all(".")

set.seed(23)
ncores <- detectCores()
cl <- makeCluster(getOption("cl.cores", ncores - 2))
doParallel::registerDoParallel(cl)

pf <- foreach(i = 1:1e3, .packages = "pfilter", .inorder = F) %dopar% {
  extract_trace(params = theta,
                y = y20210101,
                data = dat,
                data_type = "confirmation",
                npart = 1e4,
                tend = nrow(dat),
                dt = 0.2,
                error_pdf = "negbin",
                negbin_size = 20)
}
parallel::stopCluster(cl)
saveRDS(pf, "daily_sim/pf.rds")
Rt <- as.data.frame(sapply(pf, function(x) x[, "Rt"]))
daily_conf <- as.data.frame(sapply(pf, function(x) x[, "CR"]))
## save as the csv files as they are easier to
readr::write_csv(Rt, "daily_sim/Rt.csv")
readr::write_csv(daily_conf, "daily_sim/daily_confirmed.csv")
