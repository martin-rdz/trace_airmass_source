#!/bin/bash

set -x

yest=$(date +%Y%m%d --date="1 day ago")

# parse command-line arguments
while [ $# -gt 0 ]; do
    case $1 in
    "--date")
        yest=$(date --date $2 +%Y%m%d)
        shift
        shift
    ;;
    "--site")
        site=$2
        shift
        shift
    ;;
esac
done

echo $yest

source /home/traceairmass/.bashrc

downdaystart=$(date +%Y%m%d -d "${yest} -11days" )
downdayend=$(date +%Y%m%d -d "${yest} +3days" )
echo "${downdaystart} ${yest}"

echo "$USER"
cd /home/traceairmass/trace
/home/traceairmass/trace-env/bin/python3 download_gfs.py --daterange $downdaystart-$downdayend


