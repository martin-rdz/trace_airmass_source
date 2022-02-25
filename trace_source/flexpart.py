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
        for i in range(1, data['no_releases']+1):
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
        self.stat_lat = {}

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
        self.thres_categories = defaultdict(lambda: np.empty((0,)))


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



    def add_partposits_thres(self, array):
        """
        """

        for rh in self.config['height']['reception']:
            if rh == 'md':
                coords = array[array[:,3] < array[:,9]]
            else:
                coords = array[array[:,3] < float(rh)*1000]
            # print('loop trough reception heights ', rh, coords.shape)
            
            category = np.empty(coords.shape[0])
            category[coords[:,2] < -60] = 0
            category[(-60 < coords[:,2]) & (coords[:,2] < -30)] = 1
            category[(-30 < coords[:,2]) & (coords[:,2] < 0)] = 2
            category[(0 < coords[:,2]) & (coords[:,2] < 30)] = 3
            category[(30 < coords[:,2]) & (coords[:,2] < 60)] = 4
            category[60 < coords[:,2]] = 5

            category = category.astype(int)

            self.thres_categories[rh] = np.append(self.thres_categories[rh], category)

    def calc_thres_stat(self):
        """
        self.lat_names = {0: '<-60', 1: '-60..-30', 2:'-30..0', 3: '0..30', 4: '30..60', 5: '>60'}
        """

        occ_stat = namedtuple('occ_stat', 'no_below counter')

        for rh in self.config['height']['reception']:
            cat_this_height = self.thres_categories[rh]
            no = float(cat_this_height.shape[0]) if cat_this_height.shape[0] > 0 else -1
            c = {x: cat_this_height.tolist().count(x)/float(no) for x in range(6)}

            if rh != 'md':
                rh_string = rh + 'km'
            else:
                rh_string = rh

            print(rh_string, no, c)
            self.stat_lat['lat_ens_below' + rh_string] = occ_stat(no_below=no, counter=c)





class assemble_time_height(trace_source.assemble_pattern):

    #@profile
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
        days_avail = [f for f in days_avail if len(f) == 11]
        print(days_avail)
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


        self.lat_names = {0: '<-60', 1: '-60..-30', 2:'-30..0', 3: '0..30', 4: '30..60', 5: '>60'}
        self.statlat_dict = defaultdict(lambda: np.zeros((len(self.dt_list),
                                                         len(self.height_list),
                                                         len(list(self.lat_names.keys())))))


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

            self.no_part.append(traj_meta['releases_meta'][1]['no_particles'])
            self.time_res.append(10*24/len(files_for_time))

            # different structure than hysplit
            # 1. loop through the ending times of the current day
            # 2. load partposit for a specified time
            # 3. loop through heights

            for f in files_for_time:
                print('files_for_time ', f)
                part_pos = read_partpositions(folder + f, 1, ctable=True)
                part_pos = np.array(part_pos)

                for ih, h in enumerate(self.height_list):
                    #print("at ", ih, h)
                    this_population = np.where(part_pos[:,0] == ih+1)[0]
                    #release_sel = np.array([list(p) for p in part_pos if p[0]==ih+1])
                    release_sel = part_pos[this_population, :]
                    #assert np.all(release_sel == other_release)
                    meta = traj_meta['releases_meta'][ih+1]
                    #print(meta)
                    assert np.mean(meta['heights']) == h, f"{meta['heights']} {h} do not fit"
                    flex_stat[ih].add_partposits_gn(release_sel)

                    flex_stat[ih].add_partposits_ls(release_sel)
                    flex_stat[ih].add_partposits_thres(release_sel)

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

                flex_stat[ih].calc_thres_stat()
                for k in list(flex_stat[ih].stat_lat.keys()):
                    self.stat2d_dict[k+'_no_below'][it, ih] = flex_stat[ih].stat_lat[k].no_below
                    print('stat_lat ', h, k, flex_stat[ih].stat_lat[k])
                    self.statlat_dict[k][it, ih] = list(flex_stat[ih].stat_lat[k].counter.values())


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


def redistribute_rgb(r, g, b):
    threshold = 255.999
    m = max(r, g, b)
    if m <= threshold:
        return int(r), int(g), int(b)
    total = r + g + b
    if total >= 3 * threshold:
        return int(threshold), int(threshold), int(threshold)
    x = (3 * threshold - total) / (3 * m - total)
    gray = threshold - x * m
    return int(gray + x * r), int(gray + x * g), int(gray + x * b)

def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])



def plot_part_loc_map(part_pos, release_no, dt, traj, savepath, ls=None, 
        config=None, add_dyn=True, add_fire=None):
    """"""
    
    release_sel = np.array([list(p) for p in part_pos if p[0]==release_no])
    meta = traj['releases_meta'][release_no]
      
    import matplotlib
    matplotlib.use('Agg')
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    from fastkml import kml

    if ls is None:
        ls = trace_source.land_sfc.land_sfc()


    if not os.path.isdir(savepath):
        os.makedirs(savepath)
    
    if config['plotmap']['maptype'] == 'northpole':
        fig = plt.figure(figsize=(11, 10))
        if "centerlon" in config['plotmap']:
            ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=config['plotmap']['centerlon']))
        else:
            ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=45))

    if config['plotmap']['maptype'] == 'miller':
        fig = plt.figure(figsize=(10, 7))
        if "centerlon" in config['plotmap']:
            ax = plt.axes(projection=ccrs.Miller(central_longitude=config['plotmap']['centerlon']))
        else:
            ax = plt.axes(projection=ccrs.Miller(central_longitude=0))
    else:
        print("using standard map")
        fig = plt.figure(figsize=(10, 7))
        ax = plt.axes(projection=ccrs.Miller(central_longitude=0))

    assert config is not None
    if config is not None and "bounds" in config['plotmap']:
        ax.set_extent(config['plotmap']['bounds'], crs=ccrs.PlateCarree())
        
    ####
    # make a color map of fixed colors
    colors = ['lightskyblue', 'darkgreen', 'khaki', 'palegreen', 'red', 'white', 'tan']
    # better green for the light map
    colors = ['lightskyblue', 'seagreen', 'khaki', '#6edd6e', 'darkmagenta', 'white', 'tan']
    #colors = ['lightskyblue', 'seagreen', 'khaki', '#a4dd6e', 'red', 'white', 'tan']
    colors = ['lightskyblue', '#22a361', 'khaki', '#72d472', 'darkmagenta', 'white', 'tan']
    cmap = matplotlib.colors.ListedColormap(colors)

    cs = [adjust_lightness(c, amount=1.15) for c in colors]
    print('cs ', cs)
    cmap = matplotlib.colors.ListedColormap(cs)

    
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    ####
    pcm = ax.pcolormesh(ls.longs, ls.lats, ls.land_sfc_data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
    # high resolution coastlines
    #ax.coastlines(resolution='110m')
    ax.coastlines(resolution='50m')
    # ax.coastlines(resolution='10m')

    labels = ['Water', 'Forest', 'Savanna/shrubland', 'Grass/cropland', 'Urban', 'Snow/ice', 'Barren']
    handles = []
    for c, l in zip(cs, labels):
        if 'ice' in l:
            handles.append(matplotlib.patches.Patch(
                facecolor=c, label=l, edgecolor='grey', linewidth=2))
        else:
            handles.append(matplotlib.patches.Patch(color=c, label=l))
    
    # # # The modis fire map
    if add_fire:
        from cartopy.io.shapereader import Reader
        from cartopy.feature import ShapelyFeature
        # # fname = '../data/DL_FIRE_M6_78471/fire_nrt_M6_78471.shp'
        #    fname = config_dir['partposit_dir'] + '20191130_DL_FIRE_V1_90092/fire_nrt_V1_90092.shp'
        fname = f'data/fire_data/DL_FIRE_{add_fire}/fire_nrt_{add_fire}.shp'


        shp_data = Reader(fname)
        # print(next(shp_data.records()))
        frp = np.array([p.attributes['FRP'] for p in shp_data.records()])
        # print(frp.shape, np.mean(frp), np.percentile(frp, [10,25,50,75,90]))

        points = [p for p in list(shp_data.records()) if p.attributes['FRP'] > 12]
        # # for viirs pixel with their smaller size
        # points = [p for p in list(shp_data.records()) if p.attributes['FRP'] > 4]
        # #points = sorted(points, key=lambda x: x.attributes['FRP'])

        scat = ax.scatter([p.geometry.x for p in points],
                         [p.geometry.y for p in points],
                         #transform=ccrs.Geodetic(), s=0.5, 
                         transform=ccrs.PlateCarree(), s=0.5, 
                         c='crimson')

    # optionally add the dynamics from gfs grib
    if add_dyn:
        import xarray as xr

        if dt.hour%6 == 0:
            gribfile = f"data/gfs_083.2/{dt.strftime('%Y%m%d%H')}"
            print(gribfile)
            if not os.path.isfile(gribfile):
                print('try forecast file')
                gribfile = f"data/gfs_083.2/{dt.strftime('%Y%m%d%H')}_f"

            ds_isobaric = xr.load_dataset(
                gribfile, 
                engine='cfgrib', backend_kwargs={'filter_by_keys':{'typeOfLevel':'isobaricInhPa'}, 
                                                'errors': 'ignore'})

            ds_mean_sea = xr.load_dataset(
                gribfile, 
                engine='cfgrib', backend_kwargs={'filter_by_keys':{'typeOfLevel':'meanSea'}, 
                                                'errors': 'ignore'})

            prmsl = ds_mean_sea.prmsl
            gh_500 = ds_isobaric.gh.sel(isobaricInhPa=500)


        else:
            gribfile1 = f"data/gfs_083.2/{(dt - datetime.timedelta(hours=3)).strftime('%Y%m%d%H')}"
            gribfile2 = f"data/gfs_083.2/{(dt + datetime.timedelta(hours=3)).strftime('%Y%m%d%H')}"
            if not os.path.isfile(gribfile1):
                print('try forecast file')
                gribfile1 = f"data/gfs_083.2/{(dt - datetime.timedelta(hours=3)).strftime('%Y%m%d%H')}_f"
            if not os.path.isfile(gribfile2):
                print('try forecast file')
                gribfile2 = f"data/gfs_083.2/{(dt + datetime.timedelta(hours=3)).strftime('%Y%m%d%H')}_f"
            print(gribfile1, gribfile2)

            ds_isobaric = xr.load_dataset(
                gribfile1, 
                engine='cfgrib', backend_kwargs={'filter_by_keys':{'typeOfLevel':'isobaricInhPa'}, 
                                                'errors': 'ignore'})

            ds_mean_sea = xr.load_dataset(
                gribfile1, 
                engine='cfgrib', backend_kwargs={'filter_by_keys':{'typeOfLevel':'meanSea'}, 
                                                'errors': 'ignore'})

            
            ds_isobaric2 = xr.load_dataset(
                gribfile2, 
                engine='cfgrib', backend_kwargs={'filter_by_keys':{'typeOfLevel':'isobaricInhPa'}, 
                                                'errors': 'ignore'})

            ds_mean_sea2 = xr.load_dataset(
                gribfile2, 
                engine='cfgrib', backend_kwargs={'filter_by_keys':{'typeOfLevel':'meanSea'}, 
                                                'errors': 'ignore'})

            prmsl = np.mean(np.dstack((ds_mean_sea.prmsl, ds_mean_sea2.prmsl)), axis=2)
            gh_500 = np.mean(np.dstack((ds_isobaric.gh.sel(isobaricInhPa=500), 
                                        ds_isobaric2.gh.sel(isobaricInhPa=500))), axis=2)


        levels = np.arange(930,1050,4)
        cont = ax.contour(ds_mean_sea.longitude, ds_mean_sea.latitude, 
                        prmsl/100., linewidths=0.4,
                        colors='lightcoral', transform=ccrs.PlateCarree(),
                        levels=levels)

        levels = np.arange(0,900,8)
        cont = ax.contour(ds_isobaric.longitude, ds_isobaric.latitude, 
                        gh_500/10., linewidths=0.4,
                        colors='k', transform=ccrs.PlateCarree(),
                        levels=levels)

    scat = ax.scatter(release_sel[:,1], release_sel[:,2], s=2,
                      c=release_sel[:,3]/1000., cmap='plasma',
                      vmin=0.1, vmax=6.0, zorder=5,
                      #transform=ccrs.Geodetic())
                      transform=ccrs.PlateCarree())

    ax.scatter(meta['lat_lon_bounds'][0], meta['lat_lon_bounds'][1], s=14,
                      marker='^', c='tab:red', zorder=3,
                      #transform=ccrs.Geodetic())
                      transform=ccrs.PlateCarree())
    
    cbar = fig.colorbar(scat, fraction=0.025, pad=0.01)
    
    cbar.ax.set_ylabel('Height [km]', fontweight='semibold', fontsize=12)
    cbar.ax.tick_params(axis='both', which='major', labelsize=12,
                        width=2, length=4)
    
    # add the geometry boundaries...
#     k = kml.KML()
#     geo_bounds_file = '../../trace_pub/trace/data/geo_names_arctic.kml'
#     with open(geo_bounds_file) as f:
#         k.from_string(bytes(bytearray(f.read(), encoding='utf-8')))
#     docu = list(k.features())[0]
#     polygons = {}
#     colors = [(0.65098039215686276, 0.84705882352941175, 0.32941176470588235, 1.0), 
#               (1.0, 0.85098039215686272, 0.18431372549019609, 1.0), 
#               (0.89803921568627454, 0.7686274509803922, 0.58039215686274515, 1.0),
#               (0.40000000000000002, 0.76078431372549016, 0.6470588235294118, 1.0), 
#               (0.9882352941176471, 0.55294117647058827, 0.3843137254901961, 1.0), 
#               (0.55294117647058827, 0.62745098039215685, 0.79607843137254897, 1.0),  
#               (0.70196078431372544, 0.70196078431372544, 0.70196078431372544, 1.0)]
#     for p in list(docu.features()):
#         print(p.name)
#         #print(p.geometry)
#         polygons[p.name] = p.geometry
#     for i, (name, poly) in enumerate(list(polygons.items())):
#         print(i, name)
#         #ax.add_geometries([poly], ccrs.Geodetic(), alpha=0.5, facecolor=colors[i])
#         ax.add_geometries([poly], ccrs.Geodetic(), edgecolor=colors[i], facecolor='none', lw=3)

    
    ax.gridlines(linestyle=':')
    #ax.set_extent([-100, 80, 10, 80], crs=ccrs.PlateCarree())
    ##ax.set_extent([-70, 50, 20, 55], crs=ccrs.PlateCarree())
    ##ax.set_extent([-50, 40, 20, 55], crs=ccrs.PlateCarree())
    ##ax.set_extent([20, 50, 20, 40], crs=ccrs.PlateCarree())
    ##ax.set_extent([25, 35, 30, 40], crs=ccrs.PlateCarree())
    # North Pole
    #ax.set_extent([-180, 180, 45, 90], crs=ccrs.PlateCarree())


    # if maptype == 'northpole':
    #     ax.set_extent([-179, 179, 45, 90], crs=ccrs.PlateCarree())
    # elif maptype == 'southernocean':
    #     ax.set_extent([-179, 179, -75, -10], crs=ccrs.PlateCarree())
    if add_fire and 'M6' in add_fire:
        str1_pos = (.15,0.080)
        str2_pos = (.15,0.040)
    else:
        str1_pos = (.2,0.080)
    
    ax.annotate("MODIS land cover classification [Broxton and Zeng 2014, JAMC]",
                #xy=(.15, 0.108), xycoords='figure fraction',
                xy=str1_pos, xycoords='figure fraction',
                #xy=(.2, 0.085), xycoords='figure fraction',
                horizontalalignment='left', verticalalignment='bottom',
                fontsize=11)
    if add_fire and 'M6' in add_fire:
        ax.annotate("MODIS Active Fire Product, FRP > 12 MW/pixel [DOI:10.5067/FIRMS/MODIS/MCD14DL.NRT.006]",
                    #xy=(.15, 0.040), xycoords='figure fraction',
                    xy=str2_pos, xycoords='figure fraction',
                    horizontalalignment='left', verticalalignment='bottom',
                    fontsize=11)
    #ax.annotate("VIIRS Active Fire Product, FRP > 4 MW/pixel [DOI:10.5067/FIRMS/VIIRS/VNP14IMGT.NRT.001]",
    #        xy=(.2, 0.065), xycoords='figure fraction',
    #        horizontalalignment='left', verticalalignment='bottom',
    #        fontsize=11)


    fig.legend(handles, labels, ncol=4, fontsize=10,
               bbox_to_anchor=(0.20,0.185), loc='upper left')

    ax.set_title('Flexpart particle positions {}UTC\nRelease: [{:.2f} {:.2f} {:.2f} {:.2f}] {:.0f}-{:.0f}m'.format(
        dt.strftime('%Y-%m-%d %H'),
        *meta['lat_lon_bounds'], *meta['heights']),
        fontweight='semibold', fontsize=13)
    
    #fig.patch.set_facecolor('xkcd:mint green')
    #bd = matplotlib.path.Path([[0,0],[1,0],[0.5,0.5],[0,1],[0,0]])
    #ax.set_boundary(bd, transform=ax.transAxes)

    savename = savepath + "/" + "r{:0>2}_{}_{:.0f}_trajectories_map.png".format(
        release_no, dt.strftime("%Y%m%d_%H"), np.mean(meta['heights']))
    print(savename)
    fig.savefig(savename, dpi=180, transparent=False)
    plt.close('all')





if __name__ == '__main__':

    config_file = 'config_limassol.toml'
    config_file = 'config_ps122.toml'
    #dt = datetime.datetime.strptime(args.date, '%Y%m%d')
    dt = datetime.datetime(2017,9,14)
    dt = datetime.datetime(2019,10,5)
    dt_range = (dt, dt + datetime.timedelta(hours=23))
    #ath = assemble_time_height(config_file=config)
    #ath.assemble(dt_range=dt_range)
    #ath.dump2netcdf(model_str='flex')


    with open(config_file) as f:
        config = toml.loads(f.read())

    end = datetime.datetime(2019, 10, 9, 3)
    end = datetime.datetime(2019, 11, 30, 0)
    savepath = '{}{}_maps'.format(config['plot_dir'], end.strftime('%Y%m%d_%H'))

    print(savepath)
    
    #folder = '../../trace_pub/trace/flexpart_partposit/limassol/20170914_15/'
    folder = config['partposit_dir'] + '{}/'.format(end.strftime('%Y%m%d_%H'))
    dt_range = [end-datetime.timedelta(days=10), end]

    files = os.listdir(folder)
    files = sorted([f for f in files if 'partposit' in f])
    ls = trace_source.land_sfc.land_sfc()
    print(files)
    for f in files:

        for i in range(10,11):
            dt = datetime.datetime.strptime(f[10:], '%Y%m%d%H%M%S')
            part_pos = read_partpositions(folder + f, 1, ctable=False)

            traj = read_flexpart_traj_meta(folder + "trajectories.txt")
            plot_part_loc_map(part_pos, i+1, dt, traj, savepath, ls=ls, config=config)
        gc.collect()
    # convert -scale 70% -coalesce -layers Optimize -delay 20 -loop 0 `ls r11*.png | sort -r` r11.gif
    # 
