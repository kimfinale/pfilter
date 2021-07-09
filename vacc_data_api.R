library(XML)
library(RCurl)
library(tidyverse)
Sys.setlocale("LC_ALL", "Korean")
options(encoding = 'UTF-8')
# all
url <- "https://nip.kdca.go.kr/irgd/cov19stats.do?list=all"
# sido
url <- "http://nip.kdca.go.kr/irgd/cov19stats.do?list=sido"

xml_doc <- xmlTreeParse(url, useInternalNodes = TRUE, encoding = 'UTF-8')
xmldat <- getURL(url, encoding = 'UTF-8')
# call open api
xmldat <- readLines(url, encoding = 'UTF-8')
xml_doc <- xmlTreeParse(xmldat, useInternalNodes = TRUE, encoding = 'UTF-8')
xml_doc <- xmlInternalTreeParse(url, useInternalNodes = TRUE, encoding = "utf-8")

# xmldat <- getURL(url, encoding = 'UTF-8')
doc <- xmlParse (xmldat)


# [1] "<?xml version='1.0' encoding='UTF-8' ?>"
# [2] ""
# [3] ""
# [4] ""
# [5] ""
# [6] "<response>"
# [7] "<body>"
# [8] "<dataTime>2021.04.13 24:00:00</dataTime>"
# [9] "<items>"
# [10] "<item>"
# [11] "<tpcd>당일실적(A)</tpcd>"
# [12] "<firstCnt>43389</firstCnt>"
# [13] "<secondCnt>3</secondCnt>"
# [14] "</item>"
# [15] "<item>"
# [16] "<tpcd>전일누적(B)</tpcd>"
# [17] "<firstCnt>1195676</firstCnt>"
# [18] "<secondCnt>60564</secondCnt>"
# [19] "</item>"
# [20] "<item>"
# [21] "<tpcd>전체건수(C): (A)+(B)</tpcd>"
# [22] "<firstCnt>1239065</firstCnt>"
# [23] "<secondCnt>60567</secondCnt>"
# [24] "</item>"
# [25] "</items>"
# [26] "</body>"
# [27] "</response>"
# [28] ""
# [29] ""
