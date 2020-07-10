#!/usr/bin/env python3

#
# radenz@tropos.de
# 2019-09-23 inital file on lidardaten
# 2020-07-07 modified to work with the docker setup
#


import os
import datetime
import argparse
import subprocess
import string
import toml

def gen_dt_list(begin, end, step):
    print('daterange ', begin, 'to', end)
    dt_list = []
    current = begin
    while current <= end:
        dt_list.append(current)
        current += datetime.timedelta(hours=step)
    return dt_list


def download_command(template, dt):
    
    target = dt.strftime("%Y%m%d%H")

    if dt.hour % 6 == 0:
    	fcthour = '00'
    # for the non 6h interval timesteps no analysis is available
    # => use the 3h forecast
    elif dt.hour % 6 != 0 and dt.hour % 3 == 0:
        fcthour = '03'
        dt = dt - datetime.timedelta(hours=3)

    t = string.Template(WGET_TEMPLATE)
    command = t.substitute(year = dt.strftime('%Y'), \
                           month = dt.strftime('%m'), \
                           date = dt.strftime('%Y%m%d'), \
                           time = dt.strftime('%H'),\
                           fcthour = fcthour,\
                           target = target)
    return command




parser = argparse.ArgumentParser(usage="usage: %prog [date-interval]")
#parser.add_option("interval", help="interval %Y%m%d-%Y%m%d")
parser.add_argument("-m", "--model",
                  default='gfs1',
                  help="select the model to download {gfs1, gfs0p25}")
parser.add_argument('--daterange', help='date range in the format YYYYMMDD-YYYMMDD')
args = parser.parse_args()

print('args', args)

if args.model == 'gfs1':
    data_interval = 6
elif args.model == 'gfs0p25':
    data_interval = 6
else:
    raise ValueError

dt_interval = [datetime.datetime.strptime(s, "%Y%m%d") for s in args.daterange.split('-')]

dt_list = gen_dt_list(dt_interval[0], dt_interval[1], data_interval)

print(dt_list)

config = toml.load('server_logins.toml')

get_cookie_template = "wget --no-check-certificate --save-cookies auth.rda_ucar_edu --post-data='email=${mail}&passwd=${password}&action=login' https://rda.ucar.edu/cgi-bin/login"
t = string.Template(get_cookie_template)
get_cookie = t.substitute(mail= config['flexpart']['login'], password=config['flexpart']['password'])
print(get_cookie)

if args.model == 'gfs1':
    data_dir = "data/gfs_083.2/"
    WGET_TEMPLATE = '''wget -N --no-check-certificate --load-cookies auth.rda_ucar_edu -O ${target} https://rda.ucar.edu/data/ds083.2/grib2/${year}/${year}.${month}/fnl_${date}_${time}_00.grib2'''
elif args.model == 'gfs0p25':
    data_dir = "data/gfs_083.3/"
    WGET_TEMPLATE = '''wget -N --no-check-certificate --load-cookies auth.rda_ucar_edu -O ${target} https://rda.ucar.edu/data/ds083.3/${year}/${year}${month}/gdas1.fnl0p25.${date}${time}.f${fcthour}.grib2'''
else:
    raise ValueError

os.chdir(data_dir)

print('changed dir ', os.getcwd())

avail_dates = os.listdir('.')
#print "avail_dates ", avail_dates

needed_dates = list(set([dt.strftime("%Y%m%d%H") for dt in dt_list]) - set(avail_dates))
print('needed dates', needed_dates)

process = subprocess.run(get_cookie, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)
print(process.stdout)

for d in needed_dates:
    dt = datetime.datetime.strptime(d, "%Y%m%d%H")
    command = download_command(WGET_TEMPLATE, dt)

    print('downloading ', d, " using: ")
    print(command)
    process = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)
    print(process.stdout)




avail_dates = os.listdir('.')
needed_dates = list(set([dt.strftime("%Y%m%d%H") for dt in dt_list]) - set(avail_dates))
assert len(needed_dates) == 0, "not all dates that are required could be downloaded"
