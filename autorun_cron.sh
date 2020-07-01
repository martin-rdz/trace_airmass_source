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



if [ "$stations" == "all" ] || [ "$stations" == "punta" ]; then
  echo "start punta ${daterange} ${stations}" 
  #python3 test_data_avail.py --station punta --daterange $daterange
  python3 run_assemble_hysplit.py --station punta --daterange $daterange
  python3 plot2d.py --station punta --daterange $daterange
  #python3 compress_data.py --station punta --daterange 20181210-20181221
fi


if [ "$stations" == "all" ] || [ "$stations" == "pollytau" ]; then
  python3 test_data_avail.py --station pollytau --daterange $daterange
  python3 run_assemble_hysplit.py --station pollytau --daterange $daterange
  python3 plot2d.py --station pollytau --daterange $daterange
fi
