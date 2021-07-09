library(XML)
library(RCurl)
library(tidyverse)
Sys.setlocale("LC_ALL", "Korean")
options(encoding = "UTF-8")
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

# start_dt <- "20200301"
# end_dt <- "20200331"

## open api uri creation
uri <- paste0(service_url,
              paste0("?serviceKey=", service_key),
              paste0("&pageNo=", pg),
              paste0("&numOfRows=", rows),
              paste0("&startCreateDt=", start_dt),
              paste0("&endCreateDt=", end_dt))

# call open api
xml_doc <- xmlTreeParse(uri, useInternalNodes = TRUE, encoding = "UTF-8")
# XML Root Node
root_node <- xmlRoot(xml_doc)
# 오픈API 호출 결과 데이터의 개수 획득
num_rows <- as.numeric(xpathSApply(root_node, "//numOfRows", xmlValue))
# 전체 데이터의 개수 획득
total_count <- as.numeric(xpathSApply(root_node, "//totalCount", xmlValue))
# 총 오픈API 호출 횟수 계산
loop_count <- round(total_count / num_rows, 0)
## API 호출 횟수 보정
if(loop_count * num_rows < total_count){
  loop_count <- loopCount + 1
}

total_data <- data.frame()
# 오픈 API 호출을 총 오픈API 호출 횟수만큼 반복 실행
for(i in 1:loop_count){
  # 호출 URL 생성
  uri <-  paste0(service_url,
                 paste0("?serviceKey=", service_key),
                 paste0("&pageNo=", pg),
                 paste0("&numOfRows=", rows),
                 paste0("&startCreateDt=", start_dt),
                 paste0("&endCreateDt=", end_dt))

  doc <- xmlTreeParse(uri, useInternalNodes = TRUE, encoding = "UTF-8")
  root_node <- xmlRoot(doc)
  xml_data <- xmlToDataFrame(nodes = getNodeSet(root_node, '//item'))
  total_data <- rbind(total_data, xml_data)
}

total_data <- total_data %>%
  arrange("gubunEn") %>%
  distinct()

d <- dplyr::select(total_data, createDt, defCnt, gubun, localOccCnt)
d$createDt <- as.Date(d$createDt)
d$defCnt <- as.numeric(d$defCnt)
d$localOccCnt <- as.numeric(d$localOccCnt)

if (tail(d$gubun, 1) == "합계") {
  d <- d[nrow(d), c("createDt", "defCnt", "localOccCnt")]
}

dat <- readRDS("daily_sim/dat.rds") # last
names(d) <- names(dat)

## daily_confirmed = locally transmitted cases
## cumul_confirmed = local + imported cases
dat_added <- FALSE
if (tail(dat$date, 1) + 1 == d$date) {
  dat <- rbind(dat, tail(d, 1))
  dat_added <- TRUE
}
if (dat_added) {
  saveRDS(dat, "daily_sim/dat.rds")
  ## save the file also in csv format for easier handling in the shiny app
  readr::write_csv(dat, "daily_sim/dat.csv")
}
