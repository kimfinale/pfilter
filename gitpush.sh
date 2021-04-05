#!/bin/sh
cd C:/Users/jonghoon.kim/workspace/pfilter
# tstamp=$(date "+%Y%m%d")
fn1="daily_sim/dat.csv"
fn2="daily_sim/Rt.csv"
fn3="daily_sim/daily_confirmed.csv"

# echo "${fn1}"
git add "${fn1}" "${fn2}" "${fn3}"
git commit -m "Add new simulation outputs"
git push -f origin master
