#!/bin/bash

#KEYWORD=$1
LOGPATH=$2
DELAY=$3

mkdir -p "$LOGPATH"

while :
do
  # capture only processes with keyword
  #smem -t -P "(?i)$KEYWORD" > "$LOGPATH/memory_$(date --iso-8601='seconds').log"
  #top -b -c -n1 | grep -i "command\|$KEYWORD" > "$LOGPATH/cpu_$(date --iso-8601='seconds').log"

  # capture everything, process later
  COLUMNS=1000 smem -t > "$LOGPATH/memory_$(date --iso-8601='seconds').log"
  COLUMNS=1000 top -b -c -n1 -u "$USER" > "$LOGPATH/cpu_$(date --iso-8601='seconds').log"
  sleep "$DELAY"
done
