Sys.setlocale("LC_ALL", "Korean")
options(encoding = "UTF-8")

library(tidyverse)
library(rvest)

url <-  "http://ncov.mohw.go.kr/tcmBoardView.do?contSeq=369378#"

week <- 1:8 # starting the first week of Nov
weekly_confirmed <- c(13402, 13753, 17686, 22508, 28124, 38224, 40667, 35202)
weekly_severe <- c(341, 379, 443, 543, 622, 651, 467, 152)
weekly_deaths <- c(158, 207, 283, 309, 378, 288, 160, 35)

weekly_severe / weekly_confirmed
weekly_deaths / weekly_confirmed
