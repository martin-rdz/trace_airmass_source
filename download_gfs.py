#!/usr/bin/env python3

#
# radenz@tropos.de
# 2019-09-23 inital file on lidardaten
# 2020-07-07 modified to work with the docker setup
#


import os, glob
import datetime
import argparse
import subprocess
import string
import toml


# helper function for synoptic intervals
# credit https://gist.github.com/treyhunner/6218526
def divmod(dt, delta):
   seconds = int((dt - datetime.datetime.min).total_seconds())
   print('delta_total_seconds', delta.total_seconds())
   remainder = datetime.timedelta(
       seconds=int(seconds) % delta.total_seconds(),
       microseconds=dt.microsecond)
   print(dt, seconds, remainder)
   print(dt.microsecond)
   quotient = dt - remainder
   print(quotient)
   return quotient, remainder

def gen_dt_list(begin, end, step):
    print('daterange ', begin, 'to', end)
    dt_list = []
    current = begin
    while current <= end:
        dt_list.append(current)
        current += datetime.timedelta(hours=step)
    return dt_list


#command = download_command(WGET_TEMPLATE_FCST00, dt, dt_analysis=dt_last_analysis, ending_f=True)
def download_command(template, dt, fcsthour=None, ending_f=False):
    """
    Args:
        template: string
        dt: datetime
        dt_analysis:
        ending_f: include flag for forecast/fluctuating files
    """
    

    # legacy for the 3 hour gfs 0.25 files
    #if not fcst_hour and dt.hour % 6 == 0:
    # 	fcsthour = '00'
    # for the non 6h interval timesteps no analysis is available
    # => use the 3h forecast
    #elif dt.hour % 6 != 0 and dt.hour % 3 == 0:
    #    fcsthour = '03'
    #    dt = dt - datetime.timedelta(hours=3)
    if fcsthour is None:
        fcsthour = ''
        target = dt.strftime("%Y%m%d%H")
    else:
        target = (dt + datetime.timedelta(hours=fcsthour)).strftime("%Y%m%d%H")
        fcsthour = f"{fcsthour:03.0f}"
    if ending_f:
        target += "_f"


    t = string.Template(template)
    command = t.substitute(year = dt.strftime('%Y'), \
                           month = dt.strftime('%m'), \
                           date = dt.strftime('%Y%m%d'), \
                           hour = dt.strftime('%H'),\
                           fcsthour = fcsthour,\
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
    WGET_TEMPLATE = '''wget -N --no-check-certificate --load-cookies auth.rda_ucar_edu -O ${target} https://rda.ucar.edu/data/ds083.2/grib2/${year}/${year}.${month}/fnl_${date}_${hour}_00.grib2'''
    WGET_TEMPLATE_FCST00 = '''wget -N --no-check-certificate -O ${target} https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.${date}/${hour}/atmos/gfs.t${hour}z.pgrb2.1p00.f${fcsthour}'''
elif args.model == 'gfs0p25':
    data_dir = "data/gfs_083.3/"
    WGET_TEMPLATE = '''wget -N --no-check-certificate --load-cookies auth.rda_ucar_edu -O ${target} https://rda.ucar.edu/data/ds083.3/${year}/${year}${month}/gdas1.fnl0p25.${date}${hour}.f${fcsthour}.grib2'''
else:
    raise ValueError

os.chdir(data_dir)
print('changed dir ', os.getcwd())

# remove fluctuating files
for f in glob.glob("*_f"):
    print('remove', f)
    os.remove(f)
# remove small files
for f in os.listdir('.'):
    if os.path.getsize(f) < 160e3:
        print('remove', f)
        os.remove(f)


avail_dates = os.listdir('.')
avail_dates = [s.replace("_f","") for s in os.listdir('.')]

needed_dates = list(set([dt.strftime("%Y%m%d%H") for dt in dt_list]) - set(avail_dates))
print('needed dates', needed_dates)

process = subprocess.run(get_cookie, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)
print(process.stdout)


dt_now = datetime.datetime.today()
print(type(dt_now))
print(type(dt_now - datetime.timedelta(hours=4)))
# datetime of last (likely) AVAILABLE analysis
dt_last_analysis = divmod((dt_now - datetime.timedelta(hours=4)), datetime.timedelta(hours=12))[0]
print("dt_now", dt_now)
print("dt_last_analysis",dt_last_analysis)

# not as easy as initiall thought, because the ds083.2 is only updated once daily
for d in needed_dates:
    dt = datetime.datetime.strptime(d, "%Y%m%d%H")

    print('downloading ', d, " using: ")
    # for the current issues with the ncar data source
    if dt.date() < dt_now.date() - datetime.timedelta(days=3) and dt_now.hour > 12: # the database is updated quite late
        # old data, that should be in the permanent database
        print('old data')
        command = download_command(WGET_TEMPLATE, dt)
    #
    #! add the 4 hours usual runtime of gfs
    elif dt < dt_last_analysis:
        prior_analysis = divmod(dt, datetime.timedelta(hours=12))[0]
        fcsthour = (dt - prior_analysis).seconds / 3600
        print('data before last analysis', prior_analysis, fcsthour)
        # try analysis before
        command = download_command(WGET_TEMPLATE_FCST00, prior_analysis, fcsthour=fcsthour, ending_f=True)
    elif dt == dt_last_analysis:
        print('data at last analysis')
        command = download_command(WGET_TEMPLATE_FCST00, dt_last_analysis, fcsthour=0, ending_f=True)
    else:
        print(d, 'true forecast data')
        fcsthour = (dt - dt_last_analysis).total_seconds() / 3600
        print('last analysis', dt_last_analysis, 'fcsthour', fcsthour)
        command = download_command(WGET_TEMPLATE_FCST00, dt_last_analysis, fcsthour=fcsthour, ending_f=True)
    print(command)

    process = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)
    print(process.stdout)


avail_dates = [s.replace("_f","") for s in os.listdir('.')]
needed_dates = list(set([dt.strftime("%Y%m%d%H") for dt in dt_list]) - set(avail_dates))
print(needed_dates)
assert len(needed_dates) == 0, "not all dates that are required could be downloaded"

# remove small files
for f in os.listdir('.'):
    if os.path.getsize(f) < 160e3:
        print('remove', f)
        os.remove(f)
