library(XML)
library(RCurl)
library(tidyverse)
Sys.setlocale("LC_ALL", "Korean")
options(encoding = "utf-8")

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
start_dt <- "20210412"
end_dt <- gsub("-", "", as.character(Sys.Date()))
start_dt <- end_dt

## open api uri creation
uri <- paste0(service_url,
              paste0("?serviceKey=", service_key),
              paste0("&pageNo=", pg),
              paste0("&numOfRows=", rows),
              paste0("&startCreateDt=", start_dt),
              paste0("&endCreateDt=", end_dt))

# call open api
xml_doc <- xmlTreeParse(uri, useInternalNodes = TRUE, encoding = "utf-8")

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

  # 오픈 API를 호출하여 XML 데이터 획득
  doc <- xmlTreeParse(uri, useInternalNodes = TRUE, encoding = "UTF-8")
  # XML 데이터의 Root Node에 접근
  root_node <- xmlRoot(doc)
  # item Node의 데이터 추출
  xml_data <- xmlToDataFrame(nodes = getNodeSet(root_node, '//item'))
  # 추출한 데이터를 전체 데이터를 저장한 변수에 누적 저장
  total_data <- rbind(total_data, xml_data)
}



# 데이터 정렬
total_data <- total_data %>%
  arrange("gabun") %>%
  distinct()
# 데이터 확인하기
View(total_data)

# d <- dplyr::select(total_data, createDt, defCnt, gubun, localOccCnt)
# d$createDt <- as.Date(d$createDt)
# d$defCnt <- as.numeric(d$defCnt)
# d$localOccCnt <- as.numeric(d$localOccCnt)
#
# if (tail(d$gubun, 1) == "합계") {
#   d <- d[nrow(d), c("createDt", "defCnt", "localOccCnt")]
# }
#
# dat <- readRDS("daily_sim/dat.rds") # last
# names(d) <- names(dat)

# ## daily_confirmed = locally transmitted cases
# ## cumul_confirmed = local + imported cases
# dat_added <- FALSE
# if (tail(dat$date, 1) + 1 == d$date) {
#   dat <- rbind(dat, tail(d, 1))
#   dat_added <- TRUE
# }
# if (dat_added) {
#   saveRDS(dat, "daily_sim/dat.rds")
# ## save the file also in csv format for easier handling in the shiny app
#   readr::write_csv(dat, "daily_sim/dat.csv")
# }


head(totalData, 50) # 데이터 앞 부분 확인
tail(totalData, 50) # 데이터 뒷 부분 확인

# 필요한 컬럼 추출
totalData2 <- totalData %>%
  select(accDefRate, accExamCnt, accExamCompCnt, careCnt, clearCnt, deathCnt, decideCnt, examCnt, resutlNegCnt,
         stateDt, stateTime)

# 데이터 타입 변경
str(totalData2)
totalData2$accDefRate <- as.numeric(totalData2$accDefRate)
totalData2$accExamCnt <- as.numeric(totalData2$accExamCnt)
totalData2$accExamCompCnt <- as.numeric(totalData2$accExamCompCnt)
totalData2$careCnt <- as.numeric(totalData2$careCnt)
totalData2$clearCnt <- as.numeric(totalData2$clearCnt)
totalData2$deathCnt <- as.numeric(totalData2$deathCnt)
totalData2$decideCnt <- as.numeric(totalData2$decideCnt)
totalData2$examCnt <- as.numeric(totalData2$examCnt)
totalData2$resutlNegCnt <- as.numeric(totalData2$resutlNegCnt)
totalData2$stateDt <- as.Date(totalData2$stateDt, format = '%Y%m%d')

str(totalData2) # 데이터 구조 확인


## 일단위  누적 검사수 / 누적 검사 완료 수 / 검사 진행 수 / 결과 음성수 / 치료중 환자수 / 사망자 /  격리해제 / 확진자 /누적 확진률
## 데이터 집계

daily_corona <- totalData2 %>%
  group_by(stateDt) %>%
  summarise(accExamCnt = max(accExamCnt),
            accExamCompCnt = max(accExamCompCnt),
            examCnt = max(examCnt),
            resutlNegCnt = max(resutlNegCnt),
            careCnt = max(careCnt),
            deathCnt = max(deathCnt),
            clearCnt = max(clearCnt),
            decideCnt = max(decideCnt),
            accDefRate = max(accDefRate)) %>%
  arrange(stateDt)

head(daily_corona)
tail(daily_corona)

## 일별 확진자수 / 사망자수 / 격리해제수

daily_corona <- daily_corona %>%
  mutate(day_ExamCnt = accExamCnt - lag(accExamCnt, default = first(accExamCnt)),
         day_ExamCompCnt = accExamCompCnt - lag(accExamCompCnt, default = first(accExamCompCnt)),
         day_NegCnt = resutlNegCnt - lag(resutlNegCnt, default = first(resutlNegCnt)),
         day_deathCnt = deathCnt - lag(deathCnt, default = first(deathCnt)),
         day_clearCnt = clearCnt - lag(clearCnt, default = first(clearCnt)),
         day_decideCnt = decideCnt - lag(decideCnt, default = first(decideCnt)))

View(tail(daily_corona))

## == 코로나 환자 추이 == ##

## 순 확진자수 추이

ggplot(daily_corona, aes(x = stateDt, y = careCnt, group = 1)) +
  geom_line(stat = "identity") +
  geom_point() +
  scale_x_date(date_breaks = "2 week", date_labels = "%b/%d", limits = c(as.Date('2020-03-01'), NA))


## 8월 이후 신규 확진자 수 추이

daily_corona %>%
  filter(stateDt >= '2020-08-01') %>%
  ggplot(aes(x = stateDt, y = day_decideCnt, group = 1)) +
  geom_line(stat = "identity") +
  geom_point() +
  scale_x_date(date_breaks = "1 week", date_labels = "%b/%d") +
  geom_hline(yintercept = 100, color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = as.Date('2020-10-12'), color = 'blue', alpha = 0.4) +
  annotate("rect", xmin = as.Date('2020-10-05'), xmax = as.Date('2020-10-19'), ymin = 180, ymax = 210, alpha = 0.8, fill = 'white') +
  annotate("text", x = as.Date('2020-10-12'), y = 200, label = "사회적 거리두기 \n1단계 전환", fontface = 2)
