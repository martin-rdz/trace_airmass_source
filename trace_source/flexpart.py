#! /usr/bin/env python3
# coding=utf-8
""""""
"""
Author: radenz@tropos.de


""" 

import sys, os
import re
import gc
import datetime
from collections import defaultdict, Counter, namedtuple
import numpy as np
import toml
import netCDF4
import bcolz

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../')
import trace_source


def read_flexpart_traj_meta(fname, ncluster = 5):
    """ """
    
    data = {}
    with open(fname) as f:
        l = f.readline().split()
        data['end_of_sim'] = l[0] + "_" + l[1].zfill(6)
        data['version'] = l[2]
        l = f.readline().split()
        print("second line? ", l)
        l = f.readline().split()

        data['no_releases'] = int(l[0])
        data['ncluster'] = ncluster
        data['releases_meta'] = {} 
        data['releases_traj'] = defaultdict(lambda: defaultdict(list))
        for i in range(data['no_releases']):
            l = f.readline().split()
            data['releases_meta'][i] = {}
            data['releases_meta'][i]['start_times'] = list(map(float, l[0:2]))
            data['releases_meta'][i]['lat_lon_bounds'] = list(map(float, l[2:6]))
            data['releases_meta'][i]['heights'] = list(map(float, l[6:8]))
            data['releases_meta'][i]['species'] = float(l[8])
            data['releases_meta'][i]['no_particles'] = float(l[9])
            data['releases_meta'][i]['string'] =  f.readline().strip()
            #print('releases meta', data['releases_meta'][i])
            
        for line in f:
            l = line.split()
            i = int(l.pop(0))
            
            props = ['age', 'lon', 'lat', 'height', 'mean_topo',
                     'mean_mlh', 'mean_tph', 'mean_PV', 'rms_distance',
                     'rms', 'zrms_distance', 'zrms', 'frac_ml', 'frac_lt_2pvu',
                     'frac_tp']
            for p in props:
                match = re.match('([-+]?[0-9]*\.?[0-9]*)(-[0-9]*\.?[0-9]*)', l[0])
                if match and not match.group(1) == '':
                    l[0] = match.group(1)
                    l.insert(1, match.group(2))
                elem = float(l.pop(0))
                #print(p, elem)
                data['releases_traj'][i][p].append(elem)
                
                
            # cluster are not continuous
            # fix ourselfs
            avail_clusters = list(range(ncluster))
            cluster_data = []
            for k in avail_clusters:
                cluster_props = ['lon', 'lat', 'height', 'frac', 'rms']
                cluster_data.append({})
                for cp in cluster_props:
                    match = re.match('([-+]?[0-9]*\.?[0-9]*)(-[0-9]*\.?[0-9]*)', l[0])
                    if match and not match.group(1) == '':
                        l[0] = match.group(1)
                        l.insert(1, match.group(2))
                    elem = float(l.pop(0))
                    #key = 'c{}_{}'.format(k, cp)
                    #print(key, cp, elem)
                    cluster_data[k][cp] = elem
            

                for ci, elem in enumerate(cluster_data):
                    for k, v in elem.items():
                            key = 'c{}_{}'.format(ci, k)
                            data['releases_traj'][i][key].append(v)
                      
                    
            assert len(l) == 0, "line not fully consumend"
        
        return data



def get_quantized_ctable(dtype, cparams, quantize=None, expectedlen=None):
    """Return a ctable with the quantize filter enabled for floating point cols.
    
    License
        This function is taken from the reflexible package (https://github.com/spectraphilic/reflexible/tree/master/reflexible).
        Authored by John F Burkhart <jfburkhart@gmail.com> with contributions Francesc Alted <falted@gmail.com>.
        Licensed under: 'This script follows creative commons usage.'
    """
    columns, names = [], []
    for fname, ftype in dtype.descr:
        names.append(fname)
        if 'f' in ftype:
            cparams2 = bcolz.cparams(clevel=cparams.clevel, cname=cparams.cname, quantize=quantize)
            columns.append(bcolz.zeros(0, dtype=ftype, cparams=cparams2, expectedlen=expectedlen))
        else:
            columns.append(bcolz.zeros(0, dtype=ftype, cparams=cparams, expectedlen=expectedlen))
    return bcolz.ctable(columns=columns, names=names)


def read_partpositions(filename, nspec, ctable=True, clevel=5, cname="lz4", quantize=None):
    """Read the particle positions in `filename`.

    This function strives to use as less memory as possible; for this, a
    bcolz ctable container is used for holding the data.  Besides to be compressed
    in-memory, its chunked nature makes a natural fit for data that needs to
    be appended because it does not need expensive memory resize operations.

    NOTE: This code reads directly from un UNFORMATTED SEQUENTIAL data Fortran
    file so care has been taken to skip the record length at the beginning and
    the end of every record.  See:
    http://stackoverflow.com/questions/8751185/fortran-unformatted-file-format

    Parameters
    ----------
    filename : string
        The file name of the particle raw data
    nspec : int
        number of species in particle raw data
    ctable : bool
        Return a bcolz ctable container.  If not, a numpy structured array is returned instead.
    clevel : int
        Compression level for the ctable container
    cname : string
        Codec name for the ctable container.  Can be 'blosclz', 'lz4', 'zlib' or 'zstd'.
    quantize : int
        Quantize data to improve (lossy) compression.  Data is quantized using
        np.around(scale*data)/scale, where scale is 2**bits, and bits is
        determined from the quantize value.  For example, if quantize=1, bits
        will be 4.  0 means that the quantization is disabled.

    Returns
    -------
    ctable object OR structured_numpy_array

    Returning a ctable is preferred because it is used internally so it does not require to be
    converted to other formats, so it is faster and uses less memory.

    Note: Passing a `quantize` param > 0 can increase the compression ratio of the ctable
    container, but it may also slow down the reading speed significantly.

    License
        This function is taken from the reflexible package (https://github.com/spectraphilic/reflexible/tree/master/reflexible).
        Authored by John F Burkhart <jfburkhart@gmail.com> with contributions Francesc Alted <falted@gmail.com>.
        Licensed under: 'This script follows creative commons usage.'


    """

    CHUNKSIZE = 10 * 1000
    xmass_dtype = [('xmass_%d' % (i + 1), 'f4') for i in range(nspec)]
    # note age is calculated from itramem by adding itimein
    out_fields = [
                     ('npoint', 'i4'), ('xtra1', 'f4'), ('ytra1', 'f4'), ('ztra1', 'f4'),
                     ('itramem', 'i4'), ('topo', 'f4'), ('pvi', 'f4'), ('qvi', 'f4'),
                     ('rhoi', 'f4'), ('hmixi', 'f4'), ('tri', 'f4'), ('tti', 'f4')] + xmass_dtype
    raw_fields = [('begin_recsize', 'i4')] + out_fields + [('end_recsize', 'i4')]
    raw_rectype = np.dtype(raw_fields)
    recsize = raw_rectype.itemsize

    cparams = bcolz.cparams(clevel=clevel, cname=cname)
    if quantize is not None and quantize > 0:
        out = get_quantized_ctable(raw_rectype, cparams=cparams, quantize=quantize, expectedlen=int(1e6))
    else:
        out = bcolz.zeros(0, dtype=raw_rectype, cparams=cparams, expectedlen=int(1e6))

    with open(filename, "rb", buffering=1) as f:
        # The timein value is at the beginning of the file
        reclen = np.ndarray(shape=(1,), buffer=f.read(4), dtype="i4")[0]
        assert reclen == 4
        itimein = np.ndarray(shape=(1,), buffer=f.read(4), dtype="i4")
        reclen = np.ndarray(shape=(1,), buffer=f.read(4), dtype="i4")[0]
        assert reclen == 4
        nrec = 0
        while True:
            # Try to read a complete chunk
            data = f.read(CHUNKSIZE * recsize)
            read_records = int(len(data) / recsize)  # the actual number of records read
            chunk = np.ndarray(shape=(read_records,), buffer=data, dtype=raw_rectype)
            # Add the chunk to the out array
            out.append(chunk[:read_records])
            nrec += read_records
            if read_records < CHUNKSIZE:
                # We reached the end of the file
                break

    # Truncate at the max length (last row is always a sentinel, so remove it)
    out.trim(1)
    # Remove the first and last columns
    out.delcol("begin_recsize")
    out.delcol("end_recsize")

    if ctable:
        return out
    else:
        return out[:]



class flex_statistics():
    """
    build the flexpart statisctis

    this is different to the hysplit stuff, as the data format differs
    """

    def __init__(self, config, ng=None, ls=None):
        self.statistics = {}
        self.stat_ls = {}
        self.stat_gn = {}

        self.config = config
        if ng is None:
            self.ng = trace_source.land_sfc.named_geography(self.config['geonames'])
        else:
            self.ng = ng
        if ls is None:
            self.ls = trace_source.land_sfc.land_sfc()
        else:
            self.ls = ls
 

        self.gn_categories = defaultdict(lambda: np.empty((0,)))
        self.ls_categories = defaultdict(lambda: np.empty((0,)))


    def add_partposits_ls(self, array):
        """
        """

        for rh in self.config['height']['reception']:
            if rh == 'md':
                coords = array[array[:,3] < array[:,9]]
            else:
                coords = array[array[:,3] < float(rh)*1000]
            # print('loop trough reception heights ', rh, coords.shape)

            category = self.ls.get_land_sfc(coords[:,2], coords[:,1])

            self.ls_categories[rh] = np.append(self.ls_categories[rh], category)



    def calc_ls_stat(self):
        """
        """

        occ_stat = namedtuple('occ_stat', 'no_below counter')

        for rh in self.config['height']['reception']:
            cat_this_height = self.ls_categories[rh]
            no = float(cat_this_height.shape[0]) if cat_this_height.shape[0] > 0 else -1
            c = {x: cat_this_height.tolist().count(x)/float(no) for x in list(self.ls.categories.keys())}

            if rh != 'md':
                rh_string = rh + 'km'
            else:
                rh_string = rh

            print(rh_string, no, c)
            self.stat_ls['occ_ens_below' + rh_string] = occ_stat(no_below=no, counter=c)


        
    def add_partposits_gn(self, array):
        """
        """

        for rh in self.config['height']['reception']:
            if rh == 'md':
                coords = array[array[:,3] < array[:,9]]
            else:
                coords = array[array[:,3] < float(rh)*1000]
            # print('loop trough reception heights ', rh, coords.shape)

            category = self.ng.get_geo_names(coords[:,2], coords[:,1])

            self.gn_categories[rh] = np.append(self.gn_categories[rh], category)



    def calc_gn_stat(self):
        """
        """

        occ_stat = namedtuple('occ_stat', 'no_below counter')

        for rh in self.config['height']['reception']:
            cat_this_height = self.gn_categories[rh]
            no = float(cat_this_height.shape[0]) if cat_this_height.shape[0] > 0 else -1
            c = {x: cat_this_height.tolist().count(x)/float(no) for x in list(self.ng.geo_names.keys())}

            if rh != 'md':
                rh_string = rh + 'km'
            else:
                rh_string = rh

            print(rh_string, no, c)
            self.stat_gn['region_ens_below' + rh_string] = occ_stat(no_below=no, counter=c)







class assemble_time_height(trace_source.assemble_pattern):


    def assemble(self, dt_range=None):
        """
        assemble the statistics for a range of trajectories and
        save the statistics to dicts
        
        Args:
            dt_range (list(datetime), optional): timerange for that the statistics is assembled,
                default taken from config 

        """
        if dt_range is not None:
            self.dt_list = trace_source.time_list(dt_range[0],
                                                  dt_range[1],
                                                  self.config['time']['step'])

        # only for the testcase
        traj_dir = self.config['partposit_dir']
        days_avail = os.listdir(traj_dir)
        # filter only for the trajectory files with tdump extension
        folders = [f for f in days_avail if datetime.datetime.strptime(f, "%Y%m%d_%H") in self.dt_list]

        assert len(folders) > 0, 'no folders with flexpart partposit data'

        # the defaultdict is used here to sort the files by datetime within a dictionary
        # filtered_files = defaultdict(list)
        # for f in files:
        #     # regex the yyyymmdd-hh timestamp in the filename
        #     dt = datetime.datetime.strptime(re.search('([0-9]{8})-([0-9]){2}', f).group(0), '%Y%m%d-%H')
        #     height = float(re.search('([0-9]{3,6})(?=_0[0-9-]{1,4}.tdump)', f).group(0))
        #     #print(f, dt, height)
        #     if dt >= self.dt_list[0] and dt <= self.dt_list[-1]:
        #         filtered_files[dt].append((f,height))

        # here an empty dict is generated with a zero containing array
        self.stat2d_dict = defaultdict(lambda: np.zeros((len(self.dt_list), len(self.height_list))))

        self.statls_dict = defaultdict(lambda: np.zeros((len(self.dt_list), len(self.height_list), 7)))

        self.raw_dict = defaultdict(lambda: np.zeros((len(self.dt_list), len(self.height_list),
                                                     abs(self.config['time']['tr_duration'])+1)))

        # TODO make more than 7 geo names possible
        ng = trace_source.land_sfc.named_geography(self.config['geonames'])
        self.geo_names = ng.geo_names
        no_geo_names = len(list(self.geo_names.keys()))
        self.statgn_dict = defaultdict(lambda: np.zeros((len(self.dt_list),
                                                         len(self.height_list),
                                                         no_geo_names)))

        ls = trace_source.land_sfc.land_sfc()
        self.ls_categories = ls.categories


        for it, dt in enumerate(self.dt_list[:]):
            print('trajectories eding at ', dt)
            files_for_time = os.listdir(traj_dir + dt.strftime("%Y%m%d_%H"))
            files_for_time = sorted([f for f in files_for_time if "partposit_" in f])
            folder = traj_dir + dt.strftime("%Y%m%d_%H") + "/"
            print('files_for_time ', files_for_time)

            print('heights ', len(self.height_list), self.height_list)

            flex_stat = [flex_statistics(self.config, ls=ls, ng=ng) for h in self.height_list]
            traj_meta = read_flexpart_traj_meta(folder + "trajectories.txt")

            # different structure than hysplit
            # 1. loop through the ending times of the current day
            # 2. load partposit for a specified time
            # 3. loop through heights

            for f in files_for_time:
                print(f)
                part_pos = read_partpositions(folder + f, 1, ctable=False)

                for ih, h in enumerate(self.height_list):
                    print("at ", ih, h)
                    release_sel = np.array([list(p) for p in part_pos if p[0]==ih+1])
                    meta = traj_meta['releases_meta'][ih]
                    print(meta)
                    flex_stat[ih].add_partposits_gn(release_sel)

                    flex_stat[ih].add_partposits_ls(release_sel)

            # now assemble the statistics for all heights
            for ih, h in enumerate(self.height_list): 
                flex_stat[ih].calc_gn_stat()
                for k in list(flex_stat[ih].stat_gn.keys()):
                    self.stat2d_dict[k+'_no_below'][it, ih] = flex_stat[ih].stat_gn[k].no_below
                    print('stat gn ', h, k, flex_stat[ih].stat_gn[k])
                    self.statgn_dict[k][it, ih] = list(flex_stat[ih].stat_gn[k].counter.values())

                flex_stat[ih].calc_ls_stat()
                for k in list(flex_stat[ih].stat_ls.keys()):
                    self.stat2d_dict[k+'_no_below'][it, ih] = flex_stat[ih].stat_ls[k].no_below
                    print('stat ls ', h, k, flex_stat[ih].stat_ls[k])
                    self.statls_dict[k][it, ih] = list(flex_stat[ih].stat_ls[k].counter.values())


            # #assert len(f_list) > 1
            # for ih, f in enumerate(f_list):
            #     print(it, ih, f[1], dt)
            #     traj = trajectory(self.config)
            #     traj.load_file(traj_dir+f[0], silent=True)
            #     savepath = '{}/{}'.format(self.config['plot_dir'], dt.strftime('%Y%m%d'))


            #     if "timeinterval" in self.config['plotmap']:
            #         timeinterval = self.config['plotmap']['timeinterval']
            #     else:
            #         timeinterval = 12
            #     if "heights" in self.config['plotmap']:
            #         heightlist = self.config['plotmap']['heights']
            #     else:
            #         heightlist = [1500.0, 3000.0, 4500.0]
            #     #if f[1] == 3000.0 and dt.hour % 12 == 0:
            #     if f[1] in heightlist and dt.hour % timeinterval == 0:
            #         print("plotting ", f[1], dt.hour)
            #         plot_trajectories_ens(traj, savepath, ls=ls, config=self.config)
            #     #continue

            #     traj.evaluate(silent=True)
            #     traj.add_land_sfc(ls, silent=True)
            #     traj.add_ensemble_land_sfc(ls)
            #     traj.add_ensemble_geo_names(ng)
            #     #traj.add_area_land_sfc('md', ls, silent=True)
            #     #traj.add_area_land_sfc(2000, ls, silent=True)

            #     #print("at step", it, dt, ih, f)
            #     #print('keys ', traj.statistics.keys())
            #     # now the empty dict is filled with the keys (and values) of the statistics dict from traj
            #     for k in list(traj.statistics.keys()):
            #         self.stat2d_dict[k][it, ih] = traj.statistics[k]
            #     # subset of trajectory data to collect
            #     param_collect = ['latitude', 'longitude', 'height', "PRESSURE", "AIR_TEMP",
            #                      "RAINFALL", "RELHUMID", "TERR_MSL", 'age']
            #     if 'land_sfc_category' in list(traj.data.keys()):
            #         param_collect.append('land_sfc_category')
            #     for k in param_collect:
            #         #self.raw_dict[k][it, ih, :traj.data[1][k].shape[0]] = traj.data[1][k]
            #         self.raw_dict[k][it, ih, :] = traj.data[1][k]
            #         #self.raw_dict[k][it, ih, traj.data[1][k].shape[0]:] = -999.

            #     for k in list(traj.stat_ls.keys()):
            #         self.stat2d_dict[k+'_no_below'][it, ih] = traj.stat_ls[k].no_below
            #         print('stat ls ', k, traj.stat_ls[k])
            #         self.statls_dict[k][it, ih] = list(traj.stat_ls[k].counter.values())

            #     for k in list(traj.stat_gn.keys()):
            #         self.stat2d_dict[k+'_no_below'][it, ih] = traj.stat_gn[k].no_below
            #         print('stat gn ', k, traj.stat_gn[k])
            #         self.statgn_dict[k][it, ih] = list(traj.stat_gn[k].counter.values())

        # trying to free memory
        del ls
        del ng


if __name__ == '__main__':

    config = 'config_limassol.toml'
    config = 'config_ps122.toml'
    #dt = datetime.datetime.strptime(args.date, '%Y%m%d')
    dt = datetime.datetime(2017,9,14)
    dt = datetime.datetime(2019,10,5)
    dt_range = (dt, dt + datetime.timedelta(hours=23))
    ath = assemble_time_height(config_file=config)
    ath.assemble(dt_range=dt_range)
    ath.dump2netcdf(model_str='flex')
