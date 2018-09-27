#! /usr/bin/env python3
# coding=utf-8
"""
Author: radenz@tropos.de

"""

import datetime
from collections import defaultdict
import os, re
#import numpy as np
import argparse
import trace_source
import zipfile
import toml


def compress_day(filtered_files, dt_to_compress, compressdir):
    #print(filtered_files[dt_to_compress])
    
    files = filtered_files[dt_to_compress]
    assert len(files) > 0, 'no files available for this day' 

    zipfname = compressdir + '/{}_traj.zip'.format(dt_to_compress.strftime('%Y%m%d'))
    with zipfile.ZipFile(zipfname, 'w', zipfile.ZIP_DEFLATED) as myzip:
        for elem in files:
            myzip.write(elem)
    print('compression finished')

    for elem in files:
        os.remove(elem)
    print('removed original files')

parser = argparse.ArgumentParser()
parser.add_argument('--station', help='station name like limassol, barbados or mcmurdo')
parser.add_argument('--date', help='date in the format YYYYMMDD')
parser.add_argument('--daterange', help='date range in the format YYYYMMDD-YYYMMDD')

args = parser.parse_args()


if args.station is not None:
	config_file = 'config_{}.toml'.format(args.station)
else:
	config_file = 'config.toml'

with open(config_file) as config_file:
    config = toml.loads(config_file.read())
trajdir = config['traj_dir']
retdir = os.getcwd()
os.chdir(trajdir)

files = os.listdir('.')
files = [f for f in files if f[-6:] == '.tdump']
filtered_files = defaultdict(list)

for f in files[:]:
    dt = datetime.datetime.strptime(re.search('([0-9]{8})', f).group(0), '%Y%m%d')
    filtered_files[dt].append(f)
    #print(f, dt)

if args.date is not None:
    print('date to compress ', args.date)
    dt_to_compress = datetime.datetime.strptime(args.date, '%Y%m%d')
    print(dt_to_compress)
    compress_day(filtered_files, dt_to_compress, '.')

if args.daterange is not None:
    begin = datetime.datetime.strptime(args.daterange.split('-')[0], '%Y%m%d')
    end = datetime.datetime.strptime(args.daterange.split('-')[1], '%Y%m%d')
    print('compress from ', begin, 'to', end)
    current = begin
    while current <= end:
        if current in filtered_files.keys():
            print(current)
            compress_day(filtered_files, current, '.')
        else:
            print("no files for date", current)

        current += datetime.timedelta(days=1)

os.chdir(retdir)
