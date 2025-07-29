#!/usr/bin/env python3

#
# radenz@tropos.de
# 2019-09-23 inital file on lidardaten
# 2020-07-07 modified to work with the docker setup
# 2024-10-04 branced off to work with the sea ice data
#


import os, glob
import re
import datetime
import argparse
import subprocess
import string
import toml
import requests
import shutil
import traceback


import locale
locale.setlocale(locale.LC_ALL, 'en_US')



def gen_dt_list(begin, end, step):
    print('daterange ', begin, 'to', end)
    dt_list = []
    current = begin
    while current <= end:
        dt_list.append(current)
        current += datetime.timedelta(hours=step)
    return dt_list


#command = download_command(WGET_TEMPLATE_FCST00, dt, dt_analysis=dt_last_analysis, ending_f=True)
def download_url(template, dt):
    """
    Args:
        template: string
        dt: datetime
    """

    t = string.Template(template)
    url = t.substitute(year = dt.strftime('%Y'), \
                       month_str = dt.strftime('%b').lower(), \
                       date = dt.strftime('%Y%m%d'), \
                       )
    return url


parser = argparse.ArgumentParser(usage="usage: %prog [date-interval]")
#parser.add_option("interval", help="interval %Y%m%d-%Y%m%d")
parser.add_argument("-m", "--model",
                  default='gfs1',
                  help="select the model to download {gfs1, gfs0p25}")
parser.add_argument('--daterange', help='date range in the format YYYYMMDD-YYYMMDD')
args = parser.parse_args()

print('args', args)


dt_interval = [datetime.datetime.strptime(s, "%Y%m%d") for s in args.daterange.split('-')]

dt_list = gen_dt_list(dt_interval[0], dt_interval[1], 24)

print('dt_list ', dt_list)


data_dir = "data/sea_ice_amsr"
URL = '''https://data.seaice.uni-bremen.de/amsr2/asi_daygrid_swath/s6250/${year}/${month_str}/Antarctic/asi-AMSR2-s6250-${date}-v5.4.tif'''
# alternate server for the June 2025 maintenance period
URL = '''https://seaice.uni-bremen.de/data/amsr2/asi_daygrid_swath/s6250/${year}/${month_str}/Antarctic/asi-AMSR2-s6250-${date}-v5.4.tif'''

os.chdir(data_dir)
print('changed dir ', os.getcwd())

avail_dates = [re.findall(r"\d{8}",f)[0] for f in os.listdir('.')]

print('avail dates ', len(avail_dates), avail_dates[-10:])


needed_dates = list(set([dt.strftime("%Y%m%d") for dt in dt_list]) - set(avail_dates))
print('needed dates', needed_dates)


dt_now = datetime.datetime.today()
# datetime of last (likely) AVAILABLE analysis
print("dt_now", dt_now)

# not as easy as initiall thought, because the ds083.2 is only updated once daily
for d in needed_dates:
    dt = datetime.datetime.strptime(d, "%Y%m%d")

    print('downloading ', d, " using: ")
    # for the current issues with the ncar data source
    print('dt', dt)
    command = download_url(URL, dt)
    target = command.split("/")[-1]
    print(command)

    try:
        #with requests.get(command, cookies=aut.cookies, allow_redirects=True, stream=True) as r:
        with requests.get(command, allow_redirects=True, stream=True) as r:
            r.raise_for_status()
            print(r.headers)
            if 'Content-Length' in r.headers.keys():
                print('resp status', r.status_code, ' content length ', int(r.headers['Content-Length']))
                print(target)
                if int(r.headers['Content-Length']) > 10:
                    with open(target, 'wb') as f:
                        shutil.copyfileobj(r.raw, f)
            else:
                print('conent length too short')
                #with open(target, 'wb') as f:
                #    shutil.copyfileobj(r.raw, f)
    except Exception:
        print(traceback.format_exc())

avail_dates = [re.findall(r"\d{8}",f)[0] for f in os.listdir('.')]
needed_dates = list(set([dt.strftime("%Y%m%d") for dt in dt_list]) - set(avail_dates))
print(needed_dates)
assert len(needed_dates) == 0, "not all dates that are required could be downloaded"

