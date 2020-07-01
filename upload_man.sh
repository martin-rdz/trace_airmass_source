#!/bin/bash


if [ "$1" == "" ]; then
  yest=$(date +%Y%m%d --date="1 day ago")
else
  yest=$1
fi
echo $yest

year=${yest:0:4}
scp plots/punta/${yest}/{*map.png,*prof.png,*-abs-*.png} cloudnet@rsd.tropos.de:data/punta-arenas/products/trace-backtrajectories/${year}/;
scp output/punta/${yest}*.nc cloudnet@rsd.tropos.de:data/punta-arenas/products/trace-backtrajectories/${year}/;


#year=${yest:0:4}
#mon=${yest:4:2}
#day=${yest:6:2}
#ssh Picasso@rsd.tropos.de "mkdir -p pictures/trajectory_results/tel-aviv/${year}/${mon}/${day}/"
#scp plots/tel-aviv/${yest}/{*map.png,*prof.png,*-abs-*.png} Picasso@rsd.tropos.de:pictures/trajectory_results/tel-aviv/${year}/${mon}/${day}/;
##scp output/punta/${yest}*.nc cloudnet@rsd.tropos.de:data/punta-arenas/products/trace-backtrajectories/${year}/;
