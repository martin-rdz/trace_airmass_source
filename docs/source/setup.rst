####################################
Setup
####################################

site config
-----------
Config files are formatted with the `toml mini language <https://github.com/toml-lang/toml>`__.


.. code-block::

    output_dir = "output/limassol/"
    plot_dir = "plots/limassol"
    traj_dir = "trajectories/limassol/"
    partposit_dir = "flexpart_partposit/limassol/"
    
    # geonames that should be used
    # name is defined in geonames_config.toml 
    # and links the names to the .kml files
    geonames = 'standard'
    
    [time]
        # begin and end definitions only relevant if 
        # gen_hysplit_input.py is used 
        #begin = '2018-03-14_00'
        #end = '2018-03-25_00'
        begin = '2015-09-09_00'
        end = '2015-09-14_00'
        # time step for which trajectories are calculated
        step = 3
        # duration of each trajectory
        tr_duration = -240
    
    [flexpart]
        no_particles = 500
    
    [plotmap]
        maptype = 'miller'
        timeinterval = 6
        heights = [1500.0, 3000.0, 4500.0, 6000.0]
        bounds = [-120, 80, 0, 85]
        centerlon = -10
    
    [height]
        top = 12000
        plottop = 10250
        # reception heights that are used
        # 'md' mixing depth, '2.0' in km without unit
        reception = ['md', '2.0', '5.0']
    
    [station]
        name = 'Limassol, Cyprus'
        short_name = 'limassol'
        lat = 34.677
        lon = 33.038
        altitude = 10
    


FLEXPART simulations
---------------------

FLEXPART simulations are controlled with with template files in :file:`flexpart_simulations/`. Some settings are predefined in the :file:`*_template` files, `${var}` are adaped for each run with the python template library.
Per default, the particle positions during the backward simulation are dumped every 3 hours.


HYSPLIT simulations
---------------------

HYSPLIT backward trajectories are required, a 10-day 27-member ensemble setup is recommended. Please note that HYSPLIT itself is not provided within this package (The binary to run is `hyts_ens` with the respective `CONTROL` file).
Meteorological input data for HYSPLIT are taken from the GDAS1 dataset (<https://www.ready.noaa.gov/gdas1.php>) provided by the Air Resources Laboratory (ARL) of the U.S. National Weather Service’s National Centers for Environmental Prediction (NCEP).
An trajectory is calculated every 3h in steps of 500m. Conveniently the input trajectories are placed in the :file:`trajectories` directory.

Following setting provide a starting point:

.. code-block::

    SETUP.CFG
    
        &SETUP
        KMSL=0,
        tm_rain=1,
        tm_tpot=0,
        tm_tamb=1,
        tm_mixd=1,
        tm_relh=1,
        tm_terr=1,
        dxf=0.4,
        dyf=0.4,
        dzf=0.008,
        /
    
.. code-block::

    TRAJ.CFG
    
        &SETUP
        tratio = 0.75,
        delt = 0.0,
        mgmin = 10,
        khmax = 9999,
        kmixd = 0,
        kmsl = 0,
        k10m = 1,
        nstr = 0,
        mhrs = 9999,
        nver = 0,
        tout = 60,
        tm_pres = 1,
        tm_tpot = 0,
        tm_tamb = 1,
        tm_rain = 1,
        tm_mixd = 1,
        tm_relh = 1,
        tm_sphu = 0,
        tm_mixr = 0,
        tm_dswf = 0,
        tm_terr = 1,
        dxf = 0.40,
        dyf = 0.40,
        dzf = 0.01,
        messg = 'MESSAGE',
        /

