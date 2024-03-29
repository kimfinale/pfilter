rm(list=ls())
# Sys.setlocale("LC_ALL", "Korean")
d <- readr::read_csv("https://raw.githubusercontent.com/jooeungen/coronaboard_kr/master/kr_regional_daily.csv")
## get the number of locally transmitted cases (excluding "검역")
library(dplyr)
d %>% filter(region != "검역") %>%
  group_by(as.factor(date)) %>%
  summarize(confirmed = sum(confirmed)) -> d2
names(d2) <- c("date", "cumul_confirmed")
d2$date <- as.Date(as.character(d2$date), format = "%Y%m%d")
d2 %>% filter(date >= as.Date("2020-12-31")) -> d3
dat <- d3 %>% filter(date >= as.Date("2021-01-01"))
dat$daily_confirmed <- diff(d3$cumul_confirmed)
saveRDS(dat, "daily_sim/dat.rds")
## save the file also in csv format for easier handling in the shiny app
readr::write_csv(dat, "daily_sim/dat.csv")

# ## Fit starting from January 1, 2021
library(pfilter)
devtools::load_all(".")
library(parallel)
library(doParallel)
library(foreach)

y20210101 <- readRDS("outputs/y20210101.rds")

theta["betavol"] <- 0.40
usethis::use_data(theta, overwrite = T)
devtools::load_all(".")

set.seed(2)
ncores <- detectCores()
cl <- makeCluster(getOption("cl.cores", ncores - 4))
doParallel::registerDoParallel(cl)

tic <- Sys.time()
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
Sys.time() - tic
saveRDS(pf, "daily_sim/pf.rds")
Rt <- as.data.frame(sapply(pf, function(x) x[, "Rt"]))
daily_conf <- as.data.frame(sapply(pf, function(x) x[, "CR"]))
## save as the csv files as they are easier to
readr::write_csv(Rt, "daily_sim/Rt.csv")
readr::write_csv(daily_conf, "daily_sim/daily_confirmed.csv")

library(RPushbullet)
tstamp <- format(Sys.time(), "%Y-%m-%d")
pbSetup(apikey = "o.KasCi9GJjdN4Q5301FFeem2lKGdVXfCu",
        conffile = ".rpushbullet.json",
        defdev = "Samsung SM-G977N")

pbPost("note", title = paste0(tstamp, " daily pf completed!"))
