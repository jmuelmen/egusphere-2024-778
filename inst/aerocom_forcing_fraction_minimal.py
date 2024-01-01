'''Builds on the 2020 forcing decomposition paper, creating a mimimal version
without all the extra diagnostics used for testing'''
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date
import traceback
import glob
import misc
import misc.stats
import misc.plotting
import misc.geo
import misc.fileops
from misc import astro
import datetime
import cartopy
import cartopy.crs as ccrs
from shapely.geometry import MultiLineString
from copy import copy
import pickle
import progressbar

# Careful!!!
import warnings
warnings.filterwarnings("ignore")


def custom_setup():
    # Reduce the resolution of the cartopy map to make smaller size images
    feature = cartopy.feature.NaturalEarthFeature(
        'physical', 'coastline', '110m')
    f = []
    for x in list(feature.geometries()):
        b = x.bounds
        if misc.geo.haversine(b[1], b[0], b[3], b[2]) > 450:
            f.append(MultiLineString(
                [x.simplify(0.75, preserve_topology=False)]))

    cartopy.feature._NATURAL_EARTH_GEOM_CACHE[
        ('coastline', 'physical', 'custom')] = f


custom_setup()


def get_rsdt_time(utc, lon=None, lat=None, atm_corr=False):
    if lon is None:
        lat = np.arange(-88.75, 90, 2.5)
        lon = np.arange(1.25, 360, 2.5)

    try:
        doy = utc.timetuple().tm_yday
    except:
        doy = (utc -
               netcdftime._netcdftime.Datetime360Day(
                   utc.year, 1, 1)).days + 1
    utc_hour = utc.hour + utc.minute/24
        
    rsdt = np.fromfunction(
        lambda y, x: astro.insolation(
            np.radians(lat)[y.astype('int')],
            doy,
            lon[x.astype('int')]/15 + utc_hour,
            atm_corr),
        (len(lat), len(lon)))

    rsdt[rsdt<134] = np.nan
    return rsdt


class ModelData():
    def __init__(self, model, toffset, doffset, full_analysis=True, ens='r1i1p1', pd=False):
        # This sets up the arrays for storing model output. To save on memory, timesteps
        # are processed individually, requiring these arrays to accumulate output
        self.model = model
        # Time offset for each 3 hour period - required for calcualting rsdt
        # Not all models have the timestep at the time they claim
        self.toffset = toffset
        self.doffset = doffset  # Day offset (some models have a different calendar)
        self.pd = pd
        # Only for CMIP5 models (some are missing required output)
        self.full_analysis = full_analysis
        self.variables = ['cdnc', 'tcc', 'lwp', 'icc', 'lcc', 'iwp',
                          'rsut', 'rsutcs', 'rsutnoa', 'rsutcsnoa',
                          'rlut', 'rlutcs', 'od550']

        modelfolder = ('/net/seldon/disk1/Data/AEROCOM/INDIRECT3-native/' +
                       model + '/renamed/')
        self.files = {}
        for name in self.variables:
            if pd is True:
                pattern = (
                    f'{modelfolder}/' +
                    f'{name}_{model[:-5]}.nc')
                #print(pattern)
                filename = [f for f in glob.glob(pattern)]
            else:
                pattern = (
                    f'{modelfolder}/' +
                    f'{name}_PI_{model[:-5]}.nc')
                filename = [f for f in glob.glob(pattern)]
            if (len(filename) > 1) or (len(filename) == 0):
                print(pattern, filename)
                raise(IOError)
            self.files[name] = Dataset(filename[0], 'r')

        # A lot of this time calculation is required to make sure the PD and PI data match
        # and for calculating the incoming solar radiation
        Vtime = self.files[self.variables[0]].variables['time']
        vtimes = Vtime[:]
        vtimes[vtimes > 80000] = 6000
        dates = num2date(vtimes, units=Vtime.units, calendar=Vtime.calendar)
        self.doys = np.array(list(
            map(lambda x: int(x.strftime("%j")), dates)))
        self.doys[vtimes > 80000] = -1
        self.hours = np.array(list(
            map(lambda x: float(x.strftime("%H")), dates)))
        self.lon = self.files[self.variables[0]].variables['lon'][:]
        self.lat = self.files[self.variables[0]].variables['lat'][:]
        
        self.shape = self.files[self.variables[0]].variables[self.variables[0]].shape[1:]
        # Accumulated variables
        # - final char
        #  - l - low/liquid
        #  - h - high/ice
        #  - i - high/ice but not thin (IWP>thresh). This is required 
        self.outputvars = [
            # Incoming solar - can be calculated if required (1)
            'rsdt',
            # Various albedos are needed for the different components (8)
            'alb', 'albnoa',
            'albcs', 'albcsnoa',
            'albcld', 'albcldnoa',
            'albcldcs', 'albcldcsnoa',
            # These are for the LW decomposition (4)
            'rlut', 'rlutcs', 'rlutcld', 'rlutcldcs',
            # Required for forcing calculation (2)
            'TCC', 'LWP',
            # Needed for threshold calculation ('i' variables) (1)
            'IWP',
            # # Probably no longer needed - cfl, no ice gridboxes (1)
            'cfl_noi',
            # For the LWP regression against albcldl (6)
            'LWPl', 'lsLWPl', 'lsalbcldlnoa',
            'lsLWPlsq', 'lsalbcldlnoasq', 
            'lsLWPlalbcldlnoa', 
            # Required for liquid-ice decomposition (3)
            'cfl', 'cfh', 'cfi',
            # Required for liquid-ice decomposition (6)
            'albcldl', 'albcldlnoa',
            'albcldh', 'albcldhnoa',
            'albcldi', 'albcldinoa',
            # cloud-clear sky albedo - used for cloud fraction compoennt (6)
            'albcldcsl', 'albcldcslnoa',
            'albcldcsh', 'albcldcshnoa',
            'albcldcsi', 'albcldcsinoa',
            # Liquid-ice decomposition for longwave (6)
            'rlutcldl', 'rlutcldh', 'rlutcldi',
            'rlutcldcsl', 'rlutcldcsh', 'rlutcldcsi',
        ]

        self.output = {}
        for opname in self.outputvars:
            self.output[opname] = np.zeros(self.shape)
            self.output[opname+'_num'] = np.zeros(self.shape)

        self.file_index = 0        
        self.rsdt = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.__del__()
        
    def __del__(self):
        for name in self.files.keys():
            self.files[name].close()

    def read_data(self, index):
        # Read in the required variables - not all are required for
        # the forcing calculation though, Nd and AOD are only used for
        # subsequent analysis
        outputdata = {}
        for name in self.variables:
            outputdata[name] = self.files[name].variables[name][index]
            outputdata[name][np.where((outputdata[name] < -1) |
                                      (outputdata[name] > 1e15))] = np.nan

        if ('CND' in self.model) or ('anthsca' in self.model) or ('HadGEM' in self.model):
            outputdata['rsutnoa'] = outputdata['rsut']
            outputdata['rsutcsnoa'] = outputdata['rsutcs']

        if ('HadGEM' not in self.model) and ('UKESM' not in self.model):
            outputdata['cdnc'] /= 1000000
            
        # SPRINTARS appears to not have CF weighted cloud properties
        if ('SPRINTARS' not in self.model) and ('AACR' not in self.model):
            outputdata['cdnc'] /= outputdata['tcc']
            outputdata['lwp'] /= outputdata['tcc']
            outputdata['iwp'] /= outputdata['tcc']
            outputdata['lwp'][outputdata['tcc'] < 0.05] = 0
            outputdata['lwp'][outputdata['lwp'] <= 0] = 0

        datatime = num2date(
            self.files['rsut'].variables['time'][index],
            self.files['rsut'].variables['time'].units)
        outputdata['rsdt'] = get_rsdt_time(
            datatime+datetime.timedelta(
                hours=self.toffset(index)+24*self.doffset),
            self.lon,
            self.lat)

        mask = (outputdata['rsut'] == 0)
        for name in outputdata.keys():
            outputdata[name][mask] = np.nan

        outputdata['rsdtnan'] = copy(outputdata['rsdt'])
        outputdata['rsdt'] = misc.stats.zero_nans(outputdata['rsdt'])
        return outputdata

    def process_index(self, index, icclim, cflim, iwplim):
        # This function processes the data for each timestep, storing the results in
        # output arrays, if running the decomposition online, I would suggest this is the
        # part that could be implemented in the model itself (reducing required output)
        data = self.read_data(index)

        # Determine where data is liquid, ice or thick ice (iwp lim)
        liqmask = data['icc'] < icclim
        icemask = data['lcc'] < icclim
        iwpmask = (data['iwp'] > iwplim)
        clrmask = data['tcc'] < cflim
        
        procdata = {}
        # Variables copied
        procdata['cfl'] = data['lcc']
        procdata['cfh'] = data['icc']
        procdata['TCC'] = data['tcc']
        # Thick ice cloud fraction
        procdata['cfi'] = np.where(iwpmask, procdata['cfh'], 0)
        # Liquid cloud fraction is cases with no ice. Required for comparison with satellite values
        procdata['cfl_noi'] = np.where(liqmask, procdata['cfl'], np.nan)
        
        # Set values in clear-sky regions to zero
        procdata['LWP'] = data['lwp']
        procdata['LWP'][clrmask==1] = 0
        procdata['IWP'] = data['iwp']
        procdata['IWP'][clrmask==1] = 0

        # Shortwave variables
        procdata['rsdt'] = data['rsdt']
        procdata['alb'] = data['rsut']/data['rsdtnan']
        procdata['albnoa'] = data['rsutnoa']/data['rsdtnan']
        procdata['albcs'] = data['rsutcs']/data['rsdtnan']
        procdata['albcsnoa'] = data['rsutcsnoa']/data['rsdtnan']

        # Aim to catch errors in rsdt calculation (if calculated offline)
        flag = 100*((procdata['alb']>=1).sum() + (procdata['alb']<=0.05).sum())/np.isfinite(procdata['alb']).sum()
        if flag > 1:
            raise ValueError('Error in rsdt calculation')

        # Calculate cloud albedos
        procdata['albcld'] = ((procdata['alb']-procdata['albcs']*(1-procdata['TCC'])) /
                              (procdata['TCC']))
        procdata['albcld'][(procdata['albcld'] > 1) *
                           np.isfinite(procdata['albcld'])] = 1
        procdata['albcld'][clrmask==1] = np.nan

        procdata['albcldnoa'] = ((procdata['albnoa']-procdata['albcsnoa']*(1-procdata['TCC'])) /
                                 (procdata['TCC']))
        procdata['albcldnoa'][(procdata['albcldnoa'] > 1) *
                              np.isfinite(procdata['albcldnoa'])] = 1
        procdata['albcldnoa'][clrmask==1] = np.nan

        procdata['albcldcs'] = procdata['albcld']-procdata['albcs']
        procdata['albcldcsnoa'] = procdata['albcldnoa']-procdata['albcsnoa']

        # Longwave variables
        procdata['rlut'] = data['rlut']
        procdata['rlutcs'] = data['rlutcs']
        procdata['rlutcld'] = ((data['rlut']-data['rlutcs']*(1-data['tcc'])) /
                               (data['tcc']))
        procdata['rlutcld'][clrmask==1] = np.nan
        procdata['rlutcldcs'] = procdata['rlutcld'] - procdata['rlutcs']

        # LWP and CDNC for liquid cloud only
        procdata['LWPl'] = np.where(liqmask, procdata['LWP'], np.nan)

        # Removing the zero LWP points has a significant effect on the slope
        #  If LWP zeros are out, removing CDNC zeros does very little
        # We dont' care so much about the completely clear cases (they are rare in multi-year means)
        lsliqmask = liqmask*(procdata['LWP']>0)
        procdata['lsLWPl'] = np.where(lsliqmask, procdata['LWP'], np.nan)
        procdata['lsLWPlalbcldlnoa'] = np.where(lsliqmask, procdata['LWP']*procdata['albcldnoa'], np.nan)
        procdata['lsLWPlsq'] = np.where(lsliqmask, procdata['LWP']**2, np.nan)
        procdata['lsalbcldlnoasq'] = np.where(lsliqmask, procdata['albcldnoa']**2, np.nan)
        procdata['lsalbcldlnoa'] = np.where(lsliqmask, procdata['albcldnoa'], np.nan)

        # Store cloud albedos for the cloud types
        procdata['albcldl'] = np.where(liqmask, procdata['albcld'], np.nan)
        procdata['albcldh'] = np.where(icemask, procdata['albcld'], np.nan)
        procdata['albcldi'] = np.where(iwpmask, procdata['albcld'], np.nan)
        procdata['albcldcsl'] = np.where(liqmask, procdata['albcldcs'], np.nan)
        procdata['albcldcsh'] = np.where(icemask, procdata['albcldcs'], np.nan)
        procdata['albcldcsi'] = np.where(iwpmask, procdata['albcldcs'], np.nan)
        procdata['albcldlnoa'] = np.where(liqmask, procdata['albcldnoa'], np.nan)
        procdata['albcldhnoa'] = np.where(icemask, procdata['albcldnoa'], np.nan)
        procdata['albcldinoa'] = np.where(iwpmask, procdata['albcldnoa'], np.nan)
        procdata['albcldcslnoa'] = np.where(liqmask, procdata['albcldcsnoa'], np.nan)
        procdata['albcldcshnoa'] = np.where(icemask, procdata['albcldcsnoa'], np.nan)
        procdata['albcldcsinoa'] = np.where(iwpmask, procdata['albcldcsnoa'], np.nan)

        procdata['rlutcldl'] = np.where(liqmask, procdata['rlutcld'], np.nan)
        procdata['rlutcldh'] = np.where(icemask, procdata['rlutcld'], np.nan)
        procdata['rlutcldi'] = np.where(iwpmask, procdata['rlutcld'], np.nan)
        procdata['rlutcldcsl'] = np.where(liqmask, procdata['rlutcldcs'], np.nan)
        procdata['rlutcldcsh'] = np.where(icemask, procdata['rlutcldcs'], np.nan)
        procdata['rlutcldcsi'] = np.where(iwpmask, procdata['rlutcldcs'], np.nan)

        return procdata

    def append_index(self, index, icclim, cflim, iwplim):
        data = self.process_index(index, icclim, cflim, iwplim)
        for name in data.keys():
            locs = np.isfinite(data[name])
            np.add.at(self.output[name], locs, data[name][locs])
            np.add.at(self.output[name+'_num'], locs, 1)

    def get_output(self, name):
        return self.output[name]/self.output[name+'_num']

    def get_albcld_lwp_sens(self, log=False):
        num = self.output['lsLWPl_num']
        if not log:
            sxx = self.output['lsLWPlsq']-(self.output['lsLWPl']**2)/num
            syy = self.output['lsalbcldlnoasq']-(self.output['lsalbcldlnoa']**2)/num
            sxy = self.output['lsLWPlalbcldlnoa']-(self.output['lsLWPl']*self.output['lsalbcldlnoa'])/num
        else:
            sxx = self.output['lslnLWPlsq']-(self.output['lslnLWPl']**2)/num
            syy = self.output['lsalbcldlnoasq']-(self.output['lsalbcldlnoa']**2)/num
            sxy = self.output['lslnLWPlalbcldlnoa']-(self.output['lslnLWPl']*self.output['lsalbcldlnoa'])/num
        op = (sxy/sxx), (np.sqrt((1./(num-2))*((syy/sxx)-((sxy/sxx)**2))))
        return op[0], op[1], np.where(np.abs(op[0])>2*np.abs(op[1]), op[0], 0)


######################################################################
# Main program
######################################################################

if __name__=='__main__':
    icclim = 0.02
    cflim = 0.01
    iwplim = 0.0087 # 8.7gm-2 ~ tau=0.4 ~ MODIS cloud mask limit

    # model name, toffset, doffset
    models = [('ECHAM6-HAM-ETHZ_IND3', lambda x: [-0.45, -1.45][x%2], 0),
              #('ECHAM6-HAM-CND_IND3', lambda x: [-0.45, -1.45][x%2], 0),
              #('ECHAM6-HAM-anthsca0_IND3', lambda x: [-0.45, -1.45][x%2], 0),
              #('ECHAM6-HAM-anthsca1_IND3', lambda x: [-0.45, -1.45][x%2], 0),
              #('ECHAM6-HAM-anthsca1.5_IND3', lambda x: [-0.45, -1.45][x%2], 0),
              #('ECHAM6-HAM-anthsca2_IND3', lambda x: [-0.45, -1.45][x%2], 0),
              #('ECHAM6-HAM-anthsca4_IND3', lambda x: [-0.45, -1.45][x%2], 0),
              #('ECHAM63-HAM23_IND3', lambda x: [-0.45, -1.45][x%2], 0),
              #('HadGEM3_IND3', lambda x: 0.64, 0),
              #('UKESM1-AMIP-UKMO_IND3', lambda x: 0.16, 0),
              #('NCAR-CAM5.3_IND3', lambda x: -1.4, 0),
              #('NCAR-CAM5.3-MG2_IND3', lambda x: -1.4, 0),
              #('NCAR-CAM5.3-CLUBB_IND3', lambda x: -1.4, 0),
              #('NCAR-CAM5.3-CLUBB-MG2_IND3', lambda x: -1.4, 0),
              #('SPRINTARS_IND3', lambda x: 0.85, -3),
              #('SPRINTARS-KK_IND3', lambda x: 0.85, -3)
              ]

    for model in models:
        print(model[0])
        modeldata = {'PD': ModelData(*model, pd=True),
                     'PI': ModelData(*model, pd=False)}

        bar = progressbar.ProgressBar(max_value=modeldata['PD'].files['rsut'].variables['rsut'].shape[0], term_width=80)
        i = 0
        missed = 0
        while True:
            try:
                if i%25 == 0:
                    bar.update(i)
                modeldata['PD'].append_index(i, icclim=icclim, cflim=cflim, iwplim=iwplim)
                modeldata['PI'].append_index(i, icclim=icclim, cflim=cflim, iwplim=iwplim)
                if i > 365*8*5:
                    break
                i += 1
            except IndexError:
                break
            except KeyboardInterrupt:
                break
            except OverflowError:
                i += 1
                missed += 1
                continue
            except ValueError:
                i += 1
                missed += 1
                continue
        print('\nMissed: ', missed)

        # Given the processed/accumulated timestep data, calculate the forcing terms
        # This part could in theory be done online, but is really a post-processing
        # step, so I recommend doing it with output data.
        output = {}
        # First, get the PD-PI diffference in values
        for name in modeldata['PD'].outputvars:
            output[name] = modeldata['PD'].get_output(name)
            output['d'+name] = modeldata['PD'].get_output(name) - modeldata['PI'].get_output(name)

        output['num_liq'] = modeldata['PD'].output['lsLWPl_num']
        output['alb_lwp_sens'], output['alb_lwp_sens_err'], albsens = modeldata['PD'].get_albcld_lwp_sens()

        retdata = {}
        # Shortwave forcing
        retdata['F_SW'] = -output['rsdt']*output['dalb']
        # Clearsky changes
        retdata['F_dalbcs'] = -output['rsdt']*output['dalbcs']*(1-output['TCC'])
        retdata['F_dalbsurf'] = -output['rsdt']*output['dalbcsnoa']*(1-output['TCC'])
        # SW RFari (clear-sky - surface)
        retdata['F_RFari'] = retdata['F_dalbcs'] - retdata['F_dalbsurf']
        # Cloud albedo changes (total)
        retdata['F_dalbc'] = -output['rsdt']*output['dalbcld']*output['TCC']
        retdata['F_dalbcnoa'] = -output['rsdt']*output['dalbcldnoa']*output['TCC']
        # SW RFari (cloudy-skies)
        retdata['F_RFaric'] = retdata['F_dalbc'] - retdata['F_dalbcnoa']

        #########################################################
        # Now let's start calculating the liquid/ice components #
        #########################################################

        # Here is the cloud albedo part #
        #################################
        retdata['F_dalbclnoa'] = -output['rsdt']*output['dalbcldlnoa']*output['cfl']
        retdata['F_dalbchnoa'] = -output['rsdt']*output['dalbcldhnoa']*output['cfh']
        retdata['F_dalbcinoa'] = -output['rsdt']*output['dalbcldinoa']*output['cfi']
        # And the relevant forcings, attributing thin high cloud changes to low level liquid
        thin_high_adj = (retdata['F_dalbchnoa']-retdata['F_dalbcinoa'])
        retdata['F_dalbcl'] = thin_high_adj + retdata['F_dalbclnoa']
        retdata['F_dalbci'] = retdata['F_dalbcinoa']

        # LWP forcing (constant Nd) - linear regression
        retdata['F_dalbc_lwp'] = -output['rsdt']*output['cfl']*albsens*output['dlsLWPl']
        retdata['F_dalbc_lwp'][output['num_liq'] < 100] = 0
        # Twomey effect (constant LWP)
        retdata['F_dalbc_nd'] = retdata['F_dalbclnoa']-retdata['F_dalbc_lwp']
        # Distribute the ice-masked cloud forcing between the LWP and Nd components
        extra_f = retdata['F_dalbc_nd']/(retdata['F_dalbc_nd']+retdata['F_dalbc_lwp'])
        extra_f = np.clip(extra_f, 0, 1)
        retdata['F_dalbc_nd'] += extra_f*thin_high_adj
        retdata['F_dalbc_lwp'] += (1-extra_f)*thin_high_adj

        # And the cloud fraction change components #
        ############################################
        retdata['F_dcfnoa'] = -output['rsdt']*output['dTCC']*output['albcldcsnoa']
        retdata['F_dcflnoa'] = -output['rsdt']*output['dcfl']*output['albcldcslnoa']
        retdata['F_dcfhnoa'] = -output['rsdt']*output['dcfh']*output['albcldcshnoa']
        retdata['F_dcfinoa'] = -output['rsdt']*output['dcfi']*output['albcldcsinoa']
        
        # The liquid cloud fraction change has to be corrected for changes in overlying ice cloud
        # Assumes the ice and liquid cloud fraction changes are uncorrelated
        output['dcfl_corr'] = output['dcfl'] + output['dcfh']*(output['cfl']/(1-output['dcfh']))
        # Forcing from liquid cloud changes (corrected)
        retdata['F_dcfl'] = -output['rsdt']*output['dcfl_corr']*output['albcldcslnoa']
        # And make the corresponging adjustment to the ice cloud fraction forcing
        retdata['F_dcfi'] = retdata['F_dcfhnoa']+(retdata['F_dcflnoa']-retdata['F_dcfl'])

        # Mask cases with little liquid cloud? Not sure if this is necessary though
        retdata['F_dalbcl'][output['cfl'] < 0.02] = 0
        retdata['F_dalbc_lwp'][output['cfl'] < 0.02] = 0
        retdata['F_dcfl'][output['cfl'] < 0.02] = 0


        #####################
        # Longwave forcings #
        #####################
        retdata['F_LW'] = -output['drlut']
        retdata['F_dlwcs'] = -output['drlutcs']*(1-output['TCC'])
        retdata['F_dlwc'] = -output['drlutcld']*output['TCC']
        retdata['F_dlwcl'] = -output['drlutcldl']*output['cfl']
        retdata['F_dlwch'] = -output['drlutcldh']*output['cfh']
        retdata['F_dlwci'] = -output['drlutcldi']*output['cfi']

        retdata['F_dlwcf'] = -output['dTCC']*output['rlutcldcs']
        retdata['F_dlwcfl'] = -output['dcfl']*output['rlutcldcsl']
        retdata['F_dlwcfl_corr'] = -output['dcfl_corr']*output['rlutcldcsl']
        retdata['F_dlwcfh'] = -output['dcfh']*output['rlutcldcsh']

        # Adjust the LW forcing for the thin overlying cirrus effect
        retdata['F_LWcfi'] = retdata['F_dlwcfh']+(retdata['F_dlwcfl']-retdata['F_dlwcfl_corr'])

        def lwav(data):
            return misc.stats.lat_weighted_av(data, [-90, 90])

        for name in retdata.keys():
            print('{: <15} {:.2f}'.format(name, lwav(retdata[name])))
        
        

        splek
        
        #Store the output
        modelfolder = ('/net/seldon/disk1/Data/AEROCOM/INDIRECT3-native/' +
                       model[0] + '/renamed/')

        pattern = (
            f'{modelfolder}/' +
            f'rsut_{model[0][:-5]}.nc')
        filename = [f for f in glob.glob(pattern)][0]

        opfilename = f'forcing_{model[0]}_inative.nc'
        ncdf = Dataset(opfilename, 'w', 'NETCDF4')

        ncdf.title = "Aerosol forcing decomposition from AeroCom IND3 data"
        ncdf.setncattr('institution', "Space and Atmospheric Physics, Imperial College London.")
        ncdf.history = "CMIP5"
        ncdf.contact = "Edward Gryspeerdt (e.gryspeerdt@imperial.ac.uk)"
        ncdf.Conventions = "CF-1.6 "

        with Dataset(filename) as ncdf_old:
            ncdf.createDimension('lat', ncdf_old.dimensions['lat'].size)
            ncdf.createDimension('lon', ncdf_old.dimensions['lon'].size)
            ncdf.createDimension('time', 1)
            ncdf.createDimension('bnds', 2)

            Vlat = ncdf.createVariable('lat', 'f', ('lat',))
            Vlat.long_name = "latitude"
            Vlat.standard_name = "latitude"
            Vlat.units = "degrees_north"
            Vlat[:] = ncdf_old.variables['lat'][:]

            Vlon = ncdf.createVariable('lon', 'f', ('lon',))
            Vlon.long_name = "longitude"
            Vlon.standard_name = "longitude"
            Vlon.units = "degrees_east"
            Vlon[:] = ncdf_old.variables['lon'][:]

            Vtime = ncdf.createVariable('time', 'f8', ('time',))
            Vtime.long_name = "time"
            Vtime.units = "days since 1900-1-1 0:0:0"
            Vtime.calendar = "standard"
            Vtime[:] = np.array([0])

        for name in output.keys():
            Var = ncdf.createVariable(name, 'f8', ('time', 'lat', 'lon'))
            Var[:] = output[name][None, :, :].astype('float')
        ncdf.close()
