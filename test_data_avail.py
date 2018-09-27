#! /usr/bin/env python3
# coding=utf-8
"""
Author: radenz@tropos.de

"""

import datetime
from collections import defaultdict
import argparse
import toml
import os, re
import numpy as np
import trace_source


parser = argparse.ArgumentParser(description='Test if all the required files are available')
parser.add_argument('--station', help='station name like limassol, barbados or mcmurdo', required=True)
parser.add_argument('--daterange', help='date range in the format YYYYMMDD-YYYMMDD', required=True)
args = parser.parse_args()
dt_begin = datetime.datetime.strptime(args.daterange.split('-')[0], '%Y%m%d')
dt_end = datetime.datetime.strptime(args.daterange.split('-')[1], '%Y%m%d') + datetime.timedelta(hours=23)
station = args.station


with open('config_{}.toml'.format(args.station)) as config_file:
    config = toml.loads(config_file.read())

trajdir = config['traj_dir']

files = os.listdir(trajdir)
files = [f for f in files if f[-6:] == '.tdump']
filtered_files = defaultdict(list)

#dt_begin = datetime.datetime(2015, 11, 25)
#dt_end = datetime.datetime(2015, 12, 1)
dt_list = trace_source.time_list(dt_begin, dt_end, 3)

for f in files[:]:
    dt = datetime.datetime.strptime(re.search('([0-9]{8})-([0-9]){2}', f).group(0), '%Y%m%d-%H')
    height = float(re.search('([0-9]{3,6})(?=_0[0-9-]{1,4}.tdump)', f).group(0))
    size = os.stat(trajdir + '/' + f).st_size
    filtered_files[dt].append((f,height,size))
    #print(f, height, size, dt)

days_missing = set(dt_list) - set(filtered_files.keys())
print('days missing completely')
[print(dt) for dt in sorted(list(days_missing))]

for dt in sorted(list(filtered_files.keys())):
    files = [elem[0] for elem in filtered_files[dt]]
    heights = [elem[1] for elem in filtered_files[dt]]
    sizes = np.array([elem[2] for elem in filtered_files[dt]])
    #print(dt, heights)
    req_heights = np.arange(500, 12001, 500).tolist()
    allsizes = np.all(sizes > 710000)
    #print(allsizes)
    missingheights = set(req_heights) - set(heights)
    #print(missingheights)
    if missingheights:
        print('at ', dt, 'missing heights ', missingheights, allsizes)
    if not allsizes:
        print('at ', dt, 'small file ', np.where(sizes < 710000))
