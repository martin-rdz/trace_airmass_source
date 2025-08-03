#!/bin/bash

set -x

failfunction()
{
    if [ "$1" != 0 ]
    then echo "One of the commands has failed!!"
	 /usr/bin/mail -s "Task failed $2" "radenz@tropos.de" <<< "message"
         exit
    fi
}

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
downdayend=$(date +%Y%m%d -d "${yest} +1days" )
echo "${downdaystart} ${yest}"

echo "$USER"
cd /home/traceairmass/trace
#/home/traceairmass/trace-env/bin/python3 download_gfs.py --daterange $downdaystart-$downdayend


#docker exec -i trace_run bash -c "cd trace/ && python3 run_flexpart_trace.py --station $site --date $yest"
#docker exec -i trace_run bash -c "cd trace/ && python3 run_assemble.py --model flex --station $site --date $yest"
#docker exec -i trace_run bash -c "cd trace/ && python3 plot2d.py --model flex --station $site --date $yest"
#
# remove the -i as it has issues when accessed from cron
docker exec trace_run bash -c "cd trace/ && python3 run_flexpart_trace.py --station $site --date $yest"
failfunction "$?" "$site $yest run flex" 
docker exec trace_run bash -c "cd trace/ && python3 run_assemble.py --model flex --station $site --date $yest"
failfunction "$?" "$site $yest run assemble" 
docker exec trace_run bash -c "cd trace/ && python3 plot2d.py --model flex --station $site --date $yest"
failfunction "$?" "$site $yest run plot" 
#docker run -v `pwd`:/trace -v ~:/workingdir -v /data:/data -p 8890:8890 flex_container /bin/bash -c "cd trace/ && python3 plot2d.py --model flex --station $site --date $yest"

docker exec trace_run bash -c "cd trace/ && python3 animate_flex.py --station $site --datetime $yest-12 --levels 3 7 --dynamics true"
failfunction "$?" "$site $yest run animation" 

#
#
# and finally the upload
year=${yest:0:4}
ssh model_rsd2 "mkdir -p data/level1a/trace_airmass_source/${site}/${year}/;"
ssh model_rsd2 "mkdir -p plots/trace_airmass_source/${site}/${year}/;"
scp output/${site}/${yest}*.nc model_rsd2:data/level1a/trace_airmass_source/${site}/${year}/;
scp plots/${site}/${yest}_flex/*.png model_rsd2:plots/trace_airmass_source/${site}/${year}/;
scp plots/${site}/${yest}_12_maps/${yest}_12_*.gif model_rsd2:plots/trace_airmass_source/${site}/${year}/;


/usr/bin/mail -s "Task successful: $site $yest" "radenz@tropos.de" <<< "empty"
echo "finished"
exit 0


