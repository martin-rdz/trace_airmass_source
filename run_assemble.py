#!/usr/bin/env python3
# coding: utf-8

import numpy as np
import datetime
import argparse
import sys
sys.path.append("..")
import trace_source

parser = argparse.ArgumentParser()
parser.add_argument('--station', help='station name like limassol, barbados or mcmurdo')
parser.add_argument('--date', help='date in the format YYYYMMDD')
parser.add_argument('--daterange', help='date range in the format YYYYMMDD-YYYMMDD')
parser.add_argument('--model', help='model to plot, default hysplit', default='hysplit')

args = parser.parse_args()

if args.station is not None:
	config = 'config_{}.toml'.format(args.station)
else:
	config = 'config.toml'

if args.date is not None:
    print('date ', args.date)
    dt = datetime.datetime.strptime(args.date, '%Y%m%d')
    dt_range = (dt, dt + datetime.timedelta(hours=23))
    if args.model == 'hysplit':
        ath = trace_source.hysplit.assemble_time_height(config_file=config)
    elif args.model == 'flex':
        ath = trace_source.flexpart.assemble_time_height(config_file=config)
    ath.assemble(dt_range=dt_range)
    ath.dump2netcdf(model_str=args.model)


if args.daterange is not None:
    begin = datetime.datetime.strptime(args.daterange.split('-')[0], '%Y%m%d')
    end = datetime.datetime.strptime(args.daterange.split('-')[1], '%Y%m%d')
    print('daterange ', begin, 'to', end)
    current = begin
    while current <= end:
        dt_range = (current, current + datetime.timedelta(hours=23))
        if args.model == 'hysplit':
            ath = trace_source.hysplit.assemble_time_height(config_file=config)
        elif args.model == 'flex':
            ath = trace_source.flexpart.assemble_time_height(config_file=config)
        ath.assemble(dt_range=dt_range)
        ath.dump2netcdf(model_str=args.model)
        current += datetime.timedelta(days=1)

