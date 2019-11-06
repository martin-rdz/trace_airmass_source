
# coding: utf-8

import datetime, math
import os
import gc
import argparse
import numpy as np
#from Scientific.IO import NetCDF
import netCDF4
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ast
import toml

plt.rcParams['svg.fonttype'] = 'none'

# do we need viridis?
#print(os.listdir('output/'))
#print(os.listdir('output/barbados/'))

def gen_file_list(start, end, path, station, model):
    l = []
    current = start
    while current <= end:
        #l.append('output/{}_{}_hysplit-output.nc'.format(current.strftime('%Y%m%d'), station))
        #l.append('output/{1}/{0}_{1}_hysplit-output.nc'.format(current.strftime('%Y%m%d'), station))
        l.append(path + '/{0}_{1}_{2}-output.nc'.format(current.strftime('%Y%m%d'), station, model))
        current = current + datetime.timedelta(days=1)
    return l


def plot_source_2d(f, parameter, dt_list, dsp, config, savepath, config_dict, model):
    time_list = f.variables["timestamp"][:]
    dt_list = [datetime.datetime.fromtimestamp(time) for time in time_list]
    height_list = f.variables["range"][:]
    no_plots = len(dt_list)

    fig, axes = plt.subplots(1, no_plots, sharex=True, sharey=True, figsize=(12, 6))
    #fig, axes = plt.subplots(1, no_plots, sharex=True, sharey=True, figsize=(9, 6))
    ls_colors = ['lightskyblue', 'darkgreen', 'khaki', 'palegreen', 'red', 'gray', 'tan']
    ls_colors = ['lightskyblue', 'seagreen', 'khaki', '#6edd6e', 'darkmagenta', 'gray', 'tan']

    for it, dt in enumerate(dt_list):
        no_occ = f.variables[parameter + '_no_below'][it, :]
        no_occ = np.ma.masked_less(no_occ, 0)
        occ_height = f.variables[parameter][it, :, :]
        occ_height = np.ma.masked_less(occ_height, 0)
        if dsp == 'abs':
            occ_height = occ_height*no_occ[:, np.newaxis]
        occ_left = np.cumsum(occ_height, axis=1)

        categories = ast.literal_eval(f.variables[parameter].comment)

        l1 = axes[it].barh(height_list, occ_height[:, 0].T, 
                           align='center', height=0.35, color=ls_colors[0], edgecolor='none')
        l2 = axes[it].barh(height_list, occ_height[:, 1].T, left=occ_left[:, 0].T,
                           align='center', height=0.35, color=ls_colors[1], edgecolor='none')
        l3 = axes[it].barh(height_list, occ_height[:, 2].T, left=occ_left[:, 1].T,
                           align='center', height=0.35, color=ls_colors[2], edgecolor='none')
        l4 = axes[it].barh(height_list, occ_height[:, 3].T, left=occ_left[:, 2].T,
                           align='center', height=0.35, color=ls_colors[3], edgecolor='none')
        l5 = axes[it].barh(height_list, occ_height[:, 4].T, left=occ_left[:, 3].T,
                           align='center', height=0.35, color=ls_colors[4], edgecolor='none')
        l6 = axes[it].barh(height_list, occ_height[:, 5].T, left=occ_left[:, 4].T,
                           align='center', height=0.35, color=ls_colors[5], edgecolor='none')

        l7 = axes[it].barh(height_list, occ_height[:, 6].T, left=occ_left[:, 5].T,
                          align='center', height=0.35, color=ls_colors[6], edgecolor='none')

        #axes[it].set_ylim([0, 12])
        axes[it].set_ylim([0, config['height']['plottop']/1000.])
        axes[it].tick_params(axis='y', which='major', labelsize=13, 
                             width=1.5, length=3)
        axes[it].tick_params(axis='both', which='minor', width=1, length=2)
        axes[it].tick_params(axis='both', which='both', right=True, top=True,
                             direction='in')
        axes[it].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1.0))
        axes[it].set_xlabel(dt.strftime('%H:%M'), fontsize=13)
        

        if dsp == 'abs':
            if "2.0km" in parameter:
                xright = 5000
            else:
                xright = 7000

            if model == 'flex':
                xright = 4e5
            axes[it].set_xlim(right=xright)
            dsp_text = 'acc. residence time'
            if model != 'flex':
                axes[it].xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(3000))
            else:
                axes[it].xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
            axes[it].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        elif dsp == 'rel':
            xright = 1.
            axes[it].set_xlim(right=xright)
            dsp_text = 'rel. residence time'
        #axes[it].xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([0, xright]))
        axes[it].tick_params(axis='x', labeltop='off', labelbottom='off')


    axes[0].set_ylabel("Height [km]", fontweight='semibold', fontsize=14)
    axes[-1].tick_params(axis='x', labeltop='on', labelbottom='off', labelsize=11)
    plt.suptitle("{}   {}   {}".format(dt.strftime("%Y%m%d"), config_dict['station']["name"], parameter), 
                 fontweight='semibold', fontsize=15)
    if 'moving' in config_dict['station'].keys() and config_dict['station']['moving']:
        pass
    else:
        axes[-1].annotate("Endpoint: {:.1f} {:.1f} ".format(config_dict['station']["lat"], 
                                                            config_dict['station']["lon"]), 
                            xy=(.91, 0.96), xycoords='figure fraction',
                            horizontalalignment='center', verticalalignment='bottom',
                            fontsize=12)
    axes[-1].annotate("{} ".format(model), 
                      xy=(.91, 0.925), xycoords='figure fraction',
                      horizontalalignment='center', verticalalignment='bottom',
                      fontsize=12)
    axes[0].annotate('Time UTC', xy=(.5, .0),
                     xycoords='figure fraction',
                     horizontalalignment='center', verticalalignment='bottom',
                     fontsize=14, fontweight='semibold')
    axes[-1].annotate(dsp_text, xy=(.90, 0.86),
                      xycoords='figure fraction',
                      horizontalalignment='center', verticalalignment='bottom',
                      fontsize=12, fontweight='semibold')
    fig.legend((l1, l2, l3, l4, l5, l6, l7), list(categories.values()),
               #bbox_to_anchor=(0.85, 0.952),
               bbox_to_anchor=(0.03, 0.947), loc='upper left',
               ncol=4, fontsize=11, framealpha=0.8)

    #fig.set_tight_layout({'rect': [0, 0, 1, 1], 'pad': 0.1, 'h_pad': 1.5})
    #plt.tight_layout(w_pad=0.0002)
    plt.tight_layout(rect=[0, 0.02, 1, 0.93])
    plt.tight_layout(rect=[0, 0.02, 1, 0.88])
    fig.subplots_adjust(wspace=0)

    savename = savepath + "/" + dt.strftime("%Y%m%d") + "_{}_multi-land-use-{}-{}.png".format(short_name, dsp, parameter.replace('_', '-'))
    fig.savefig(savename, dpi=400)
    savename = savepath + "/" + dt.strftime("%Y%m%d") + "_{}_multi-land-use-{}-{}.svg".format(short_name, dsp, parameter.replace('_', '-'))
    fig.savefig(savename)
    plt.close()


def plot_geonames_2d(f, parameter, dt_list, dsp, config, savepath, config_dict, model):
    with open('geonames_config.toml') as config_file:
        geo_config = toml.loads(config_file.read())

    geo_names = {int(k): v for k, v in geo_config[config['geonames']]['geo_names'].items()}
    #assert geo_names == ast.literal_eval(f.variables[parameter].comment), "geonames not matching {}".format(str(geo_names))

    NUM_COLORS = 7
    cm = plt.get_cmap('Set2')
    cNorm  = matplotlib.colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
    colors = [scalarMap.to_rgba(i) for i in range(NUM_COLORS)]
    colors = [(0.65098039215686276, 0.84705882352941175, 0.32941176470588235, 1.0), 
              (1.0, 0.85098039215686272, 0.18431372549019609, 1.0), 
              (0.89803921568627454, 0.7686274509803922, 0.58039215686274515, 1.0),
              (0.40000000000000002, 0.76078431372549016, 0.6470588235294118, 1.0), 
              (0.9882352941176471, 0.55294117647058827, 0.3843137254901961, 1.0), 
              (0.55294117647058827, 0.62745098039215685, 0.79607843137254897, 1.0),  
              (0.70196078431372544, 0.70196078431372544, 0.70196078431372544, 1.0)]

    time_list = f.variables["timestamp"][:]
    dt_list = [datetime.datetime.fromtimestamp(time) for time in time_list]    
    height_list = f.variables["range"][:]
    no_plots = len(dt_list)

    fig, axes = plt.subplots(1, no_plots, sharex=True, sharey=True, figsize=(12, 6))
    #fig, axes = plt.subplots(1, no_plots, sharex=True, sharey=True, figsize=(9, 6))

    for it, dt in enumerate(dt_list):
        no_occ = f.variables[parameter + '_no_below'][it, :]
        no_occ = np.ma.masked_less(no_occ, 0)
        occ_height = f.variables[parameter][it, :, :]
        occ_height = np.ma.masked_less(occ_height, 0)
        if dsp == 'abs':
            occ_height = occ_height*no_occ[:, np.newaxis]
        occ_left = np.cumsum(occ_height, axis=1)

        l = []
        for i in range(len(geo_names)):
            if i == 0:
                l.append(axes[it].barh(height_list, occ_height[:, 0].T, left=0, 
                    align='center', height=0.3, color=colors[0], edgecolor='none'))
            else:
                l.append(axes[it].barh(height_list, occ_height[:, i].T, left=occ_left[:, i-1].T,
                    align='center', height=0.3, color=colors[i], edgecolor='none'))

        axes[it].set_ylim([0, 12])
        axes[it].set_ylim([0, config['height']['plottop']/1000.])
        axes[it].tick_params(axis='y', which='major', labelsize=14, 
                             width=1.5, length=3)
        axes[it].tick_params(axis='both', which='minor', width=1, length=2)
        axes[it].tick_params(axis='both', which='both', right=True, top=True,
                             direction='in')
        axes[it].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1.0))
        axes[it].set_xlabel(dt.strftime('%H:%M'), fontsize=14)
        

        if dsp == 'abs':
            if "2.0km" in parameter:
                xright = 5000
            else:
                xright = 7000
            if model == 'flex':
                xright = 4e5
            axes[it].set_xlim(right=xright)
            dsp_text = 'acc. residence time'
            if model != 'flex':
                axes[it].xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(3000))
            else:
                axes[it].xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
            axes[it].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        elif dsp == 'rel':
            xright = 1.
            axes[it].set_xlim(right=xright)
            dsp_text = 'rel. residence time'
        #axes[it].xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([0, xright]))
        axes[it].tick_params(axis='x', labeltop='off', labelbottom='off')


    axes[0].set_ylabel("Height [km]", fontweight='semibold', fontsize=14)
    axes[-1].tick_params(axis='x', labeltop='on', labelbottom='off', labelsize=11)
    plt.suptitle("{}   {}   {}".format(dt.strftime("%Y%m%d"), config_dict['station']["name"], parameter), 
                 fontweight='semibold', fontsize=15)
    if 'moving' in config_dict['station'].keys() and config_dict['station']['moving']:
        pass
    else:
        axes[-1].annotate("Endpoint: {:.1f} {:.1f} ".format(config_dict['station']["lat"], 
                                                            config_dict['station']["lon"]), 
                        xy=(.91, 0.96), xycoords='figure fraction',
                        horizontalalignment='center', verticalalignment='bottom',
                        fontsize=12)
    axes[-1].annotate("{} ".format(model), 
                      xy=(.91, 0.925), xycoords='figure fraction',
                      horizontalalignment='center', verticalalignment='bottom',
                      fontsize=12)
    axes[0].annotate('Time UTC', xy=(.5, .0),
                     xycoords='figure fraction',
                     horizontalalignment='center', verticalalignment='bottom',
                     fontsize=14, fontweight='semibold')
    axes[-1].annotate(dsp_text, xy=(.90, 0.86),
                      xycoords='figure fraction',
                      horizontalalignment='center', verticalalignment='bottom',
                      fontsize=12, fontweight='semibold')

    fig.legend(l, list(geo_names.values()),
               loc='upper left',
               #bbox_to_anchor=(0.01, 0.952),
               bbox_to_anchor=(0.01, 0.947),
               ncol=5, fontsize=12, framealpha=0.8)
    #fig.set_tight_layout({'rect': [0, 0, 1, 1], 'pad': 0.1, 'h_pad': 1.5})
    #plt.tight_layout(w_pad=0.0002)

    #plt.tight_layout(rect=[0, 0.02, 1, 0.93])
    #plt.tight_layout(rect=[0, 0.02, 1, 0.90])
    plt.tight_layout(rect=[0, 0.02, 1, 0.88])
    fig.subplots_adjust(wspace=0)

    savename = savepath + "/" + dt.strftime("%Y%m%d") + "_{}_multi-geonames-{}-{}.png".format(short_name, dsp, parameter.replace('_', '-'))
    fig.savefig(savename, dpi=400)
    savename = savepath + "/" + dt.strftime("%Y%m%d") + "_{}_multi-geonames-{}-{}.svg".format(short_name, dsp, parameter.replace('_', '-'))
    fig.savefig(savename)
    plt.close()


def plot_source_height_profile(f, parameter, dt, it, config, savepath, config_dict, model):

    dsp = 'abs'

    time_list = f.variables["timestamp"][:]
    dt_list = [datetime.datetime.fromtimestamp(time) for time in time_list]
    height_list = f.variables["range"][:]
    no_occ = f.variables[parameter + '_no_below'][it, :]
    no_occ = np.ma.masked_less(no_occ, 0)
    #print(no_occ)

    occ_height = f.variables[parameter][it, :, :]
    occ_height = np.ma.masked_less(occ_height, 0)
    #print(occ_height)

    if dsp == 'abs':
        occ_height = occ_height*no_occ[:, np.newaxis]

    occ_left = np.cumsum(occ_height, axis=1)
    #print(occ_left)

    ls_colors = ['lightskyblue', 'darkgreen', 'khaki', 'palegreen', 'red', 'grey', 'tan']
    ls_colors = ['lightskyblue', 'seagreen', 'khaki', '#6edd6e', 'darkmagenta', 'gray', 'tan']

    fig, ax = plt.subplots(1, figsize=(5, 7.5))

    ax.barh(height_list, occ_height[:, 0].T, 
            align='center', height=0.3, color=ls_colors[0], edgecolor='none')
    ax.barh(height_list, occ_height[:, 1].T, left=occ_left[:, 0].T,
            align='center', height=0.3, color=ls_colors[1], edgecolor='none')
    ax.barh(height_list, occ_height[:, 2].T, left=occ_left[:, 1].T,
            align='center', height=0.3, color=ls_colors[2], edgecolor='none')
    ax.barh(height_list, occ_height[:, 3].T, left=occ_left[:, 2].T,
            align='center', height=0.3, color=ls_colors[3], edgecolor='none')
    ax.barh(height_list, occ_height[:, 4].T, left=occ_left[:, 3].T,
            align='center', height=0.3, color=ls_colors[4], edgecolor='none')
    ax.barh(height_list, occ_height[:, 5].T, left=occ_left[:, 4].T,
            align='center', height=0.3, color=ls_colors[5], edgecolor='none')

    ax.barh(height_list, occ_height[:, 6].T, left=occ_left[:, 5].T,
            align='center', height=0.3, color=ls_colors[6], edgecolor='none')

    ax.set_ylim([0, config['height']['plottop']/1000.])
    #ax.set_ylim([0, 10.25])

    if dsp == 'abs':
        if "2.0km" in parameter:
            xright = 5000
        else:
            xright = 7000

        if model == 'flex':
            xright = 4e5
        ax.set_xlim(right=xright)
        if model != 'flex':
            ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(3000))
        else:
            ax.xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
        ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    else:
        ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.1))

    if dsp == 'abs':
        no_xloc = xright*1.13
    else:
        no_xloc = 1.17
    for i, h in enumerate(height_list):
        if h < 10:
            ax.text(no_xloc, h-0.15, int(no_occ.filled(0)[i]),
                    verticalalignment='bottom', horizontalalignment='right',
                    fontsize=11, fontweight='semibold')

    ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.5))

    ax.set_title("{} {}".format(dt.strftime("%Y%m%d_%H"), parameter), 
                 fontweight='semibold', fontsize=15)
    ax.set_xlabel("Acc. residence time", fontweight='semibold', fontsize=14)
    ax.set_ylabel("Height [km]", fontweight='semibold', fontsize=14)

    ax.tick_params(axis='both', which='major', labelsize=14, 
                   width=2, length=4)
    ax.tick_params(axis='both', which='minor', width=1.5, length=3)
    ax.tick_params(axis='both', which='both', right=True, top=True,
                   direction='in')

    plt.tight_layout(rect=[0,0,0.92,1])
    savename = savepath + "/" + dt.strftime("%Y%m%d_%H") + "_{}_land-use-{}-{}.png".format(short_name, dsp, parameter.replace('_', '-'))
    fig.savefig(savename, dpi=250)
    plt.close()


def plot_filename(filename, config_dict, model, config_file='config.toml'):
    f = netCDF4.Dataset(filename, 'r')
    time_list = f.variables["timestamp"][:]
    height_list = f.variables["range"][:]

    with open(config_file) as config_file:
        config = toml.loads(config_file.read())

    #print(time_list)
    time_delta = (time_list[1:] - time_list[:-1])[0]
    dt_list = [datetime.datetime.fromtimestamp(time) for time in time_list]
    #print(dt_list)

    #savepath = 'plots/{}/'.format(dt_list[0].strftime('%Y%m%d'))
    if model != 'hysplit':
        savepath = '{}/{}_{}/'.format(
            config['plot_dir'], dt_list[0].strftime('%Y%m%d'), model)
    else:
        savepath = '{}/{}/'.format(config['plot_dir'], dt_list[0].strftime('%Y%m%d'))
    if not os.path.isdir(savepath):
        os.makedirs(savepath)

    ## coordinate arrays for pcmesh
    #print(height_list)
    height_pcmesh = (height_list - np.mean((height_list[1:] - height_list[:-1])/2.)).tolist()
    height_pcmesh.append(height_pcmesh[-1] + np.mean((height_list[1:] - height_list[:-1])))
    #print(height_pcmesh)
    dt_pcmesh = [dt-datetime.timedelta(seconds=int(time_delta)/2) for dt in dt_list]
    dt_pcmesh.append(dt_pcmesh[-1]+datetime.timedelta(seconds=int(time_delta)))
    #print(dt_pcmesh)

    print(f)
    #print(f.variables.keys())


    if 'min_height_ag' in f.variables.keys():
        min_height_ag = f.variables['min_height_ag'][:]/1000.

        fig, ax = plt.subplots(1, figsize=(10, 5.7))
        pcmesh = ax.pcolormesh(matplotlib.dates.date2num(dt_pcmesh[:]),
                            height_pcmesh[:],
                            np.transpose(min_height_ag),
                            cmap='gist_rainbow')
        cbar = fig.colorbar(pcmesh)
        ax.set_xlim([dt_list[0], dt_list[0]+datetime.timedelta(hours=24)])
        ax.set_ylim([0, height_list[-1]])
        ax.set_xlabel("Time UTC", fontweight='semibold', fontsize=15)
        ax.set_ylabel("Height [km]", fontweight='semibold', fontsize=15)
        cbar.ax.set_ylabel("Trajectory min Height [km]", fontweight='semibold', fontsize=15)
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))

        ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=3))
        ax.xaxis.set_minor_locator(matplotlib.dates.MinuteLocator(byminute=[0,30]))
        ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(500))

        ax.tick_params(axis='both', which='major', labelsize=14, 
                    width=2, length=4)
        ax.tick_params(axis='both', which='minor', width=1.5, length=3)
        ax.tick_params(axis='both', which='both', right=True, top=True,
                    direction='in')
        cbar.ax.tick_params(axis='both', which='major', labelsize=14,
                            width=2, length=4)

        savename = savepath + "/" + dt_list[0].strftime("%Y%m%d_%H%M") + "_{}_min_height.png".format(short_name)
        fig.savefig(savename, dpi=250)
        plt.close()

    if 'Tmin_whole' in f.variables.keys():
        var = f.variables['Tmin_whole'][:]

        fig, ax = plt.subplots(1, figsize=(10, 5.7))
        pcmesh = ax.pcolormesh(matplotlib.dates.date2num(dt_pcmesh[:]),
                            height_pcmesh[:],
                            np.transpose(var),
                            vmin=-55, vmax=5,
                            cmap='gist_rainbow')
        cbar = fig.colorbar(pcmesh)
        ax.set_xlim([dt_list[0], dt_list[0]+datetime.timedelta(hours=24)])
        ax.set_ylim([0, height_list[-1]])
        ax.set_xlabel("Time UTC", fontweight='semibold', fontsize=15)
        ax.set_ylabel("Height [km]", fontweight='semibold', fontsize=15)
        cbar.ax.set_ylabel("Trajectory min T [C]", fontweight='semibold', fontsize=15)
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))

        ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=3))
        ax.xaxis.set_minor_locator(matplotlib.dates.MinuteLocator(byminute=[0,30]))
        ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(500))

        ax.tick_params(axis='both', which='major', labelsize=14, 
                    width=2, length=4)
        ax.tick_params(axis='both', which='minor', width=1.5, length=3)
        ax.tick_params(axis='both', which='both', right=True, top=True,
                    direction='in')
        cbar.ax.tick_params(axis='both', which='major', labelsize=14,
                            width=2, length=4)

        savename = savepath + "/" + dt_list[0].strftime("%Y%m%d_%H%M")+ "_{}_min_T_whole.png".format(short_name)
        fig.savefig(savename, dpi=250)
        plt.close()

    if 'Tmin_24h' in f.variables.keys():
        var = f.variables['Tmin_24h'][:]

        fig, ax = plt.subplots(1, figsize=(10, 5.7))
        pcmesh = ax.pcolormesh(matplotlib.dates.date2num(dt_pcmesh[:]),
                            height_pcmesh[:],
                            np.transpose(var),
                            vmin=-55, vmax=5,
                            cmap='gist_rainbow')
        cbar = fig.colorbar(pcmesh)
        ax.set_xlim([dt_list[0], dt_list[0]+datetime.timedelta(hours=24)])
        ax.set_ylim([0, height_list[-1]])
        ax.set_xlabel("Time UTC", fontweight='semibold', fontsize=15)
        ax.set_ylabel("Height [km]", fontweight='semibold', fontsize=15)
        cbar.ax.set_ylabel("Trajectory min T 24h [C]", fontweight='semibold', fontsize=15)
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))

        ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=3))
        ax.xaxis.set_minor_locator(matplotlib.dates.MinuteLocator(byminute=[0,30]))
        ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(500))

        ax.tick_params(axis='both', which='major', labelsize=14, 
                    width=2, length=4)
        ax.tick_params(axis='both', which='minor', width=1.5, length=3)
        ax.tick_params(axis='both', which='both', right=True, top=True,
                    direction='in')
        cbar.ax.tick_params(axis='both', which='major', labelsize=14,
                            width=2, length=4)

        savename = savepath + "/" + dt_list[0].strftime("%Y%m%d_%H%M") + "_{}_mint_T_24h.png".format(short_name)
        fig.savefig(savename, dpi=250)
        plt.close()

    if 'mean_bearing_from_endpoint' in f.variables.keys():
        var = f.variables['mean_bearing_from_endpoint'][:]

        fig, ax = plt.subplots(1, figsize=(10, 5.7))
        pcmesh = ax.pcolormesh(matplotlib.dates.date2num(dt_pcmesh[:]),
                            height_pcmesh[:],
                            np.transpose(var),
                            vmin=0, vmax=360,
                            cmap='gist_rainbow')
        cbar = fig.colorbar(pcmesh)
        ax.set_xlim([dt_list[0], dt_list[0]+datetime.timedelta(hours=24)])
        ax.set_ylim([0, height_list[-1]])
        ax.set_xlabel("Time UTC", fontweight='semibold', fontsize=15)
        ax.set_ylabel("Height [km]", fontweight='semibold', fontsize=15)
        cbar.ax.set_ylabel("Bearing of Trajectory (mean)", fontweight='semibold', fontsize=15)
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))

        ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=3))
        ax.xaxis.set_minor_locator(matplotlib.dates.MinuteLocator(byminute=[0,30]))
        ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(500))

        ax.tick_params(axis='both', which='major', labelsize=14, 
                    width=2, length=4)
        ax.tick_params(axis='both', which='minor', width=1.5, length=3)
        ax.tick_params(axis='both', which='both', right=True, top=True,
                    direction='in')
        cbar.ax.tick_params(axis='both', which='major', labelsize=14,
                            width=2, length=4)

        savename = savepath + "/" + dt_list[0].strftime("%Y%m%d_%H%M")+ "_{}_mean_bearing.png".format(short_name)
        fig.savefig(savename, dpi=250)
        plt.close()

    if 'traj_distance_total' in f.variables.keys():
        var = f.variables['traj_distance_total'][:]

        fig, ax = plt.subplots(1, figsize=(10, 5.7))
        pcmesh = ax.pcolormesh(matplotlib.dates.date2num(dt_pcmesh[:]),
                            height_pcmesh[:],
                            np.transpose(var),
                            vmin=2500, vmax=18000,
                            cmap='gist_rainbow')
        cbar = fig.colorbar(pcmesh)
        ax.set_xlim([dt_list[0], dt_list[0]+datetime.timedelta(hours=24)])
        ax.set_ylim([0, height_list[-1]])
        ax.set_xlabel("Time UTC", fontweight='semibold', fontsize=15)
        ax.set_ylabel("Height [km]", fontweight='semibold', fontsize=15)
        cbar.ax.set_ylabel("Total length [km]", fontweight='semibold', fontsize=15)
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))

        ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=3))
        ax.xaxis.set_minor_locator(matplotlib.dates.MinuteLocator(byminute=[0,30]))
        ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(500))

        ax.tick_params(axis='both', which='major', labelsize=14, 
                    width=2, length=4)
        ax.tick_params(axis='both', which='minor', width=1.5, length=3)
        ax.tick_params(axis='both', which='both', right=True, top=True,
                    direction='in')
        cbar.ax.tick_params(axis='both', which='major', labelsize=14,
                            width=2, length=4)

        savename = savepath + "/" + dt_list[0].strftime("%Y%m%d_%H%M") + "_{}_length_tot.png".format(short_name)
        fig.savefig(savename, dpi=250)
        plt.close()


    for rh in config['height']['reception']:
            if rh != 'md':
                rh_string = rh + 'km'
            else:
                rh_string = rh
            plot_source_2d(f, 'occ_ens_below'+rh_string, dt_list, 'abs', config, savepath, config_dict, model)
            plot_source_2d(f, 'occ_ens_below'+rh_string, dt_list, 'rel', config, savepath, config_dict, model)
            plot_geonames_2d(f, 'region_ens_below'+rh_string, dt_list, 'abs', config, savepath, config_dict, model)
            plot_geonames_2d(f, 'region_ens_below'+rh_string, dt_list, 'rel', config, savepath, config_dict, model)

    for dt in dt_list:
        it = dt_list.index(dt)
        #plot_source_height_profile(f, 'occ_belowmd', dt, it, config, savepath, config_dict)
        plot_source_height_profile(f, 'occ_ens_belowmd', dt, it, config, savepath, config_dict, model)
        #plot_source_height_profile(f, 'occ_below2.0km', dt, it, config, savepath, config_dict)
        plot_source_height_profile(f, 'occ_ens_below2.0km', dt, it, config, savepath, config_dict, model)

    gc.collect()

parser = argparse.ArgumentParser(description='Plot trace output')
parser.add_argument('--station', help='station name like limassol, barbados or mcmurdo', required=True)
parser.add_argument('--daterange', help='date range in the format YYYYMMDD-YYYMMDD', required=True)
parser.add_argument('--model', help='model to plot, default hysplit', default='hysplit')

args = parser.parse_args()


begin = datetime.datetime.strptime(args.daterange.split('-')[0], '%Y%m%d')
end = datetime.datetime.strptime(args.daterange.split('-')[1], '%Y%m%d')
station = args.station

if args.station is not None:
    config = 'config_{}.toml'.format(args.station)

config_dict = toml.load(config)
short_name = config_dict['station']['short_name']

#filelist = gen_file_list(datetime.date(2016, 1, 22), datetime.date(2016, 1, 24), station)
filelist = gen_file_list(begin, end, config_dict['output_dir'], short_name, args.model)
for filename in filelist:
    print('=============================================================================')
    print('plotting file ', filename)
    plot_filename(filename, config_dict, args.model, config_file=config)

