#!/bin/bash

if [ "$1" == "" ]; then
  yest=$(date +%Y%m%d --date="1 day ago")
else
  yest=$1
fi

if [ "$2" == "" ]; then
  stations="all"
else
  stations=$2
fi

echo $yest
daterange=${yest}-${yest}
echo $daterange
echo $stations

python3 run_flexpart_trace.py --station $stations --date $yest
python3 run_assemble.py --model flex --station $stations --daterange $daterange
python3 plot2d.py --model flex --station $stations --daterange $daterange

# python3 download_gfs.py --daterange 20210210-20210226
# python3 run_flexpart_trace.py --station leipzig --date 20210222
# python3 run_assemble.py --model flex --station leipzig --daterange 2021022-20210222
# python3 run_assemble.py --model flex --station leipzig --daterange 20210222-20210222
# python3 plot2d.py --model flex --station leipzig --daterange 20210222-20210222
# python3 run_flexpart_trace.py --station leipzig --date 20210221
# python3 plot2d.py --model flex --station leipzig --daterange 20210221-20210221
# python3 run_assemble.py --model flex --station leipzig --daterange 20210221-20210221
# python3 plot2d.py --model flex --station leipzig --daterange 20210221-20210221
# python3 run_flexpart_trace.py --station leipzig --date 20210222
# python3 run_flexpart_trace.py --station leipzig --date 20210220
