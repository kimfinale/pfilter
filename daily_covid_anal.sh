#!/bin/sh
cd C:/Users/jonghoon.kim/workspace/pfilter
# run the r script for download the data and run the particle filter
Rscript -e "eval(parse('daily_covid_anal.R', encoding='UTF-8'))"
# select files to push to the GitHub
fn1="daily_sim/dat.csv"
fn2="daily_sim/Rt.csv"
fn3="daily_sim/daily_confirmed.csv"
git add "${fn1}" "${fn2}" "${fn3}"
git commit -m "Add new simulation outputs"
git push -f origin master
# send the message for the simulation suing Pushbullet
Rscript -e "eval(parse('daily_covid_anal_push.R', encoding='UTF-8'))"
