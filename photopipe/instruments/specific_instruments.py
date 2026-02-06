import re
import os
from glob import glob
import sys
from datetime import datetime as dt

import numpy as np
from astropy.io import fits as pf
from astropy.units import pixel_scale
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from astropy import units as u

from .instrument_class import instrument


def calculate_airmass(hdr):
    location = EarthLocation.from_geodetic(hdr['LONGITUD'] * u.deg, hdr['LATITUDE'] * u.deg, hdr['ALTITUDE'] * u.meter)
    coords = SkyCoord(hdr['RA-D'], hdr['DEC-D'], unit=(u.deg, u.deg))
    time = (hdr['JD-BEG'] + hdr['JD-END']) / 2
    time = Time(time, format='jd', scale='utc')
    altaz = coords.transform_to(AltAz(location=location, obstime=time))
    return altaz.alt.deg, altaz.az.deg, altaz.secz.value


class ratir(instrument):
    def __init__(self):
        instrument.__init__(self, 'ratir', 4)

    def possible_filters(self):
        filters = ['u', 'g', 'r', 'i', 'Z', 'Y', 'J', 'H']
        return filters

    def has_cam_bias(self, idx):
        cam_bias = [True, True, False, False]
        return cam_bias[idx]

    def has_cam_dark(self, idx):
        cam_dark = [True, True, False, False]
        return cam_dark[idx]

    def is_cam_split(self, idx):
        CAM_SPLIT = [False, False, True, True]
        return CAM_SPLIT[idx]
        
    def change_header_keywords(self, h, cam):
        cam_i = int(cam[1])
        
        # WCS relevant parameters (RATIR H2RGs have barrel distortions)
        a = -19.60381671
        b = -4128.15179797
        CAM_SECPIX1 = [0.3168, 0.3171, 0.2988, 0.2983]
        CAM_SECPIX2 = [0.3171, 0.3191, 0.2955, -0.2945]
        CAM_THETA = [0.60, 2.40, -88.1, 90.4]
        CAM_X0 = [512, 512, 1177, 924]
        CAM_Y0 = [512, 512, 1031, 982]
        CAM_PXSCALE = [0.32, 0.32, 0.3, 0.3]  # C0, C1, C2, C3 in arcsec/px

        H2RG_ASTR = {'PV1_1': 1.0, 'PV2_1': 1.0, 'PV1_17': a, 'PV2_17': a, 'PV1_19': 2.0*a,
                     'PV2_19': 2.0*a, 'PV1_21': a, 'PV2_21': a, 'PV1_31': b, 'PV2_31': b,
                     'PV1_33': 3.0*b, 'PV2_33': 3.0*b, 'PV1_35': 3.0*b, 'PV2_35': 3.0*b,
                     'PV1_37': b, 'PV2_37': b}
                      
        # Keyword names
        RA_KEY = 'ETRRQRA'
        DEC_KEY = 'ETRRQDE'
        CENTER_KEY = 'STRRQAP'  # RATIR header keyword specifying which H2RG filters the target is focused on
        OFFRA_KEY = 'ETRRQRAO'
        OFFDEC_KEY = 'ETRRQDEO'
        SOFTGAIN_KEY = 'SOFTGAIN'
        AIRMASS_KEY = 'STROBAM'
        DATEOBS_KEY = 'SDATE'
        
        try:
            h[RA_KEY]
        except:
            print('No ' + RA_KEY + ' keyword set, use 0 for RA and offset values, set center key to rcenter')
            h[RA_KEY] = 0
            h[OFFRA_KEY] = 0
            h[OFFDEC_KEY] = 0
            h[CENTER_KEY] = 'rcenter'
        
        # Definitions for camera
        CAM_WAVE = ['OPT', 'OPT', 'IR', 'IR']
        CAM_SPLIT = [False, False, True, True]
            
        # frame corners in arcmin offset from center.  top-left, bottom-left, bottom-right, top-right.
        # top==north, left==east
        CAMOFFS = np.array([
            [[2.785, 2.632], [2.604, -2.775], [-2.800, -2.615], [-2.635, 2.789]],  # C0 corner offsets in arcmin
            [[2.817, 2.624], [2.607, -2.818], [-2.807, -2.624], [-2.616, 2.818]],   # C1 corner offsets in arcmin
            [[5.012, 6.229], [4.569, -3.905], [-0.227, -3.678], [0.228, 6.453]],    # C2-Z corner offsets in arcmin
            [[-0.556, 6.488], [-1.013, -3.642], [-5.488, -3.445], [-5.018, 6.683]],   # C2-Y corner offsets in arcmin
            [[4.834, 5.430], [4.720, -4.701], [-0.324, -4.647], [-0.185, 5.494]],   # C3-J corner offsets in arcmin
            [[-0.916, 5.503], [-1.059, -4.639], [-5.318, -4.594], [-5.154, 5.557]]])  # C3-H corner offsets in arcmin
        FRAMECENTER = CAMOFFS.mean(axis=1)  # field centers in arcmin (E,N) offset from center
        
        # Offsets of the pointing apertures east and north in arcmin
        APOFFS = {
            "rcenter":      np.array([0, 0]),
            "icenter":      np.array([0, 0]),
            "ricenter":     np.array([0, 0]),
            "riZJcenter":   np.array([1.2, 0]),
            "riYHcenter":   np.array([-1.8, 0]),
            "ZJcenter":     np.array([2.2, 0.7]),
            "YHcenter":     np.array([-3.2, 0.7])}

        try:
            prpslid = h['PRPSLID']      
            vstid = h['VSTID']
                            
            targname = '{}-vis{}'.format(prpslid, vstid) 
        except:
            sys.exit('TARGET NAME UNKNOWN Check header')
        
        if cam in ['C2a', 'C2b', 'C3a', 'C3b']:
            for key in H2RG_ASTR:
                h[key] = H2RG_ASTR[key] 
        
        SOFTGAIN = h['SOFTGAIN']
        CAM_GAIN = [16.80/SOFTGAIN, 18.64/SOFTGAIN, 2.2/SOFTGAIN, 2.4/SOFTGAIN]
                
        # set keyword values
        h['CAMERA'] = cam_i
        h['TARGNAME'] = targname
        h['OBJECT'] = targname
        h['OBJNAME'] = targname
        h['PIXSCALE'] = CAM_PXSCALE[cam_i]
        h['WAVELENG'] = CAM_WAVE[cam_i]
        h['AIRMASS'] = h[AIRMASS_KEY]
        h['DATE-OBS'] = h[DATEOBS_KEY]
        h['GAIN'] = (self.get_cam_gain(h, cam_i), 'in electrons/DN')
        h['SATURATE'] = (self.get_cam_sat(h, cam_i), 'in electrons/DN')
        h['CRPIX1'] = CAM_X0[cam_i]
        h['CRPIX2'] = CAM_Y0[cam_i]
        h['CTYPE1'] = 'RA---TAN'
        h['CTYPE2'] = 'DEC--TAN'
        h['CD1_1'] = -CAM_SECPIX1[cam_i]*np.cos(CAM_THETA[cam_i]*np.pi/180.0)/3600.
        h['CD2_1'] = CAM_SECPIX1[cam_i]*np.sin(CAM_THETA[cam_i]*np.pi/180.0)/3600.
        h['CD1_2'] = CAM_SECPIX2[cam_i]*np.sin(CAM_THETA[cam_i]*np.pi/180.0)/3600.
        h['CD2_2'] = CAM_SECPIX2[cam_i]*np.cos(CAM_THETA[cam_i]*np.pi/180.0)/3600.
        
        h['CRVAL1'] = h[RA_KEY] - APOFFS[h[CENTER_KEY]][0]/60.0 + h[OFFRA_KEY]
        # includes aperture offsets and target offsets (ie. dithering)
        h['CRVAL2'] = h[DEC_KEY] - APOFFS[h[CENTER_KEY]][1]/60.0 + h[OFFDEC_KEY]
            
        h['FILTER'] = self.get_filter(h, cam)
        h['NAXIS1'] = self.slice(cam)[1].stop - self.slice(cam)[1].start
        h['NAXIS2'] = self.slice(cam)[0].stop - self.slice(cam)[0].start
        h['CRPIX1'] = CAM_X0[cam_i] - self.slice(cam)[1].start
        h['CRPIX2'] = CAM_Y0[cam_i] - self.slice(cam)[0].start        
    
        return h
        
    def slice(self, cam):
        C0_SLICE = np.s_[0:1023, 125:1023]
        C1_SLICE = np.s_[75:1000, 15:1000]
        Z_SLICE = np.s_[4:975, 100:2000]
        Y_SLICE = np.s_[1135:2043, 240:2043]
        J_SLICE = np.s_[50:1000, 4:2000]
        H_SLICE = np.s_[1200:2043, 4:1940]
        
        slicedict = {
            'C0': C0_SLICE, 'C1': C1_SLICE, 'C2': np.s_[0:2048, 0:2048], 'C2a': Z_SLICE,
            'C2b': Y_SLICE, 'C3': np.s_[0:2048, 0:2048], 'C3a': J_SLICE, 'C3b': H_SLICE}
        
        return slicedict[cam]
        
    def get_cam_sat(self, h, idx):
        # saturation levels for each detector in DNs as a function of the SOFTGAIN keyword extracted from a frame's
        # header
        SOFTGAIN = h['SOFTGAIN']
        CAM_SAT = [(2.**16/SOFTGAIN)-1, (2.**16/SOFTGAIN)-1, (36000./SOFTGAIN)-1, (36000./SOFTGAIN)-1]

        return CAM_SAT[idx]
    
    def get_cam_gain(self, h, idx):
        SOFTGAIN = h['SOFTGAIN']
        CAM_GAIN = [16.80/SOFTGAIN, 18.64/SOFTGAIN, 2.2/SOFTGAIN, 2.4/SOFTGAIN]
        
        return CAM_GAIN[idx]        

    def get_exptime(self, h):
        return h['EXPTIME']

    def get_filter(self, h, cam):
    
        SLICE_FILTERS = {'C2': 'NA', 'C2a': 'Z', 'C2b': 'Y', 'C3': 'NA', 'C3a': 'J', 'C3b': 'H'}
        
        idx = int(cam[1])

        if self.is_cam_split(idx):
            return SLICE_FILTERS[cam]
        else:
            return h['FILTER']      
            
    def get_centered_filter(self, h, idx):

        if self.is_cam_split(idx):
            CENTER_KEY = 'STRRQAP'  # RATIR header keyword specifying which H2RG filters the target is focused on
    
            try:
                center = h[CENTER_KEY].split('center')[0]
                # RATIR header keyword specifying which H2RG filters the target is focused on
            except:
                sys.exit('CENTER keyword not found in header')  
            
            if h[CENTER_KEY].split('center')[0].count('ri') != 0:  
                center_cams = list(center)
                return center_cams[idx]
            else:
                return center

        else:
            return h['FILTER']           

    def change_file_names(self, files):
        # Has correct file format, no changes needed
        return
        
    def original_file_format(self):
        file_format = '????????T??????C??.fits'
        return file_format


class lmi(instrument):

    def __init__(self):
        instrument.__init__(self, 'lmi', 1)

    def possible_filters(self):
        filters = ['U', 'B', 'V', 'R', 'I', 'YISH', 'SDSS-U', 'SDSS-G', 'SDSS-R', 'SDSS-I', 'SDSS-Z']
        return filters

    def has_cam_bias(self, idx):
        cam_bias = [True]
        return cam_bias[idx]

    def has_cam_dark(self, idx):
        cam_dark = [False]
        return cam_dark[idx]

    def is_cam_split(self, idx):
        CAM_SPLIT = [False]
        return CAM_SPLIT[idx]
        
    def change_header_keywords(self, h, cam):
    
        f1 = h["FILTER1"]
        f2 = h["FILTER2"]
        
        if f1 == "OPEN":
            h['FILTER'] = f2
        elif f2 == "OPEN":
            h['FILTER'] = f1
        else:
            h['FILTER'] = "%s-%s" % (f1, f2)
                
        if h['OBJECT'] == '':
            h['TARGNAME'] = h['OBSTYPE']
            h['OBJECT'] = h['OBSTYPE']
            h['OBJNAME'] = h['OBSTYPE']
        else:
            h['TARGNAME'] = h['OBJNAME']
        
        # set keyword values
        h['CAMERA'] = 0
        h['BINNING'] = h['CCDSUM']
        h['BINX'] = h['ADELX_01']
        h['BINY'] = h['ADELY_01']
        h['PIXSCALE'] = 0.12*h['BINX']
        h['WAVELENG'] = 'OPT'
        h['SATURATE'] = (self.get_cam_sat(h, 'C0'), 'in electrons/DN') 
        h['UTC'] = h['UT']
        # h['NAXIS1']   = 3096-27+1
        # h['NAXIS2']   = 3079-3+1
    
        return h
        
    def slice(self, cam):
        C0_SLICE = np.s_[3:3079, 27:3096]
        return C0_SLICE
        
    def get_cam_sat(self, h, idx):
        sat = 65000.
        return sat
    
    def get_cam_gain(self, h, idx):
        gain = h['GAIN']
        
        return gain      

    def get_exptime(self, h):
        return h['EXPTIME']

    def get_filter(self, h, cam):
    
        f1 = h["FILTER1"]
        f2 = h["FILTER2"]
        
        if f1 == "OPEN":
            h['FILTER'] = f2
        elif f2 == "OPEN":
            h['FILTER'] = f1
        else:
            h['FILTER'] = "%s-%s" % (f1, f2)
            
        return h['FILTER']      
            
    def get_centered_filter(self, h, idx):

        f1 = h["FILTER1"]
        f2 = h["FILTER2"]
        
        if f1 == "OPEN":
            h['FILTER'] = f2
        elif f2 == "OPEN":
            h['FILTER'] = f1
        else:
            h['FILTER'] = "%s-%s" % (f1, f2)
            
        return h['FILTER']  

    def change_file_names(self, files):
        # Has correct _file format, no changes needed
        
        obstype_postdict = {
            'SKY FLAT': 'f', 'DOME FLAT': 'f',
            'BIAS': 'b', 'OBJECT': 'o'}
        
        datesearch = r"(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})\.\d+"

        fits_check = glob('????????T??????C??.fits')
        
        for _file in files:
            if _file in fits_check:
                continue
            print(_file)
            pyim = pf.open(_file)
            h = pyim[0].header
            if h['OBSTYPE'].upper() not in obstype_postdict.keys():
                print('Found file, {},  with OBSTYPE {} which may be an error'.format(_file, h['OBSTYPE'].upper()))
                continue
            obstype_post = obstype_postdict[h['OBSTYPE'].upper()]
            
            idate = h['DATE-OBS']
            match = re.search(datesearch, idate)
            match.groups()
            
            newname = ''
            
            for i, item in enumerate(match.groups()):
            
                if i == 3:
                    newname += 'T'
                
                newname += item
            
            newname += 'C0'+obstype_post+'.fits'
            print(newname)
            os.rename(_file, newname)
        
        return
        
    def original_file_format(self):
        file_format = 'lmi.????.fits'
        return file_format        


class rimas(instrument):

    def __init__(self):
        instrument.__init__(self, 'rimas', 2)

    def possible_filters(self):
        filters = ['Y', 'J', 'H', 'K']
        return filters

    def has_cam_bias(self, idx):
        cam_bias = [False, False]
        return cam_bias[idx]

    def has_cam_dark(self, idx):
        cam_dark = [False, False]
        return cam_dark[idx]

    def is_cam_split(self, idx):
        CAM_SPLIT = [False, False]
        return CAM_SPLIT[idx]

    def gen_wcs(self, hdr):
        image_width = hdr['NAXIS1']
        image_height = hdr['NAXIS2']
        # pixscale = hdr['PIXSCALE'] / 3600
        pixscale = hdr.get('PIXSCALE', 0.19)
        plate_scale = pixscale * u.arcsec
        # image_width = 925  # for earlier iterations such as 20251009
        # image_height = 925  # for earlier iterations such as 20251009

        w = WCS(naxis=hdr['NAXIS'])
        w.wcs.crpix = [image_width / 2, image_height / 2]
        coords = SkyCoord(hdr['RA-D'], hdr['DEC-D'], unit=(u.deg, u.deg))
        coords = coords.spherical_offsets_by(d_lon=(75 + 13) * plate_scale, d_lat=(-25 - 8) * plate_scale)  # TODO: do this offset at the telescope level
        w.wcs.crval = [coords.ra.deg, coords.dec.deg]
        # theta = np.deg2rad(5.3)  # for earlier iterations such as 20251009
        theta = np.deg2rad(hdr['ROT'])
        pixscale = pixscale / 3600
        # w.wcs.cdelt = [pixscale, pixscale]
        w.wcs.cd = [[-pixscale * np.cos(theta), pixscale * np.sin(theta)],
                    [pixscale * np.sin(theta), pixscale * np.cos(theta)]]
        w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        return w

# 653274189
    def change_header_keywords(self, h, cam):
        # set keyword values
        h['WAVELENG'] = 'IR'
        h['GAIN'] = h['GAIN0']  # TODO: fix this in the asdetector header
        h['EXPTIME'] = float(h['EXPTIME'])  # Needs to be float, can change later
        h['SATURATE'] = float(h['SATURATE'])
        h['ALT'], az, airmass = calculate_airmass(h)  # TODO: fix ALT value from labview
        h['AZ'] = h['AZIMUTH']
        h['SECZ'] =  h.get('AIRMASS', airmass)
        # idate = h['DATEOBS']  # DATEOBS needs to be DATE-OBS and use Lmi format
        # idatet = idate[0:4] + '-' + idate[4:6] + '-' + idate[6:8] + 'T'
        # h['DATE-OBS'] = idatet + idate[9:11] + ':' + idate[11:13] + ':' + idate[13:15] + '.' + idate[15:17]
        h['DATE-OBS'] = h['DATE']
        h['DATEOBS'] = h['DATE']
        h['OBSDATE'] = h['DATE']
        h['CAMERA'] = int(h['ASICINDX'])
        h['OBJECT'] = h['OBJNAME']
        h['OBSTYPE'] = h['OBJTYPE']

        h['FILTER'] = self.get_filter(h, h['CAMERA'])

        if h['OBJECT'] == '':
            h['TARGNAME'] = h['OBSTYPE']
            h['OBJECT'] = h['OBSTYPE']
            h['OBJNAME'] = h['OBSTYPE']
        else:
            h['TARGNAME'] = h['OBJNAME']
        wcs = self.gen_wcs(h)
        wcs_header = wcs.to_header(relax=True)
        wcs_header['CD1_1'] = wcs_header['PC1_1']
        wcs_header['CD2_2'] = wcs_header['PC2_2']
        del wcs_header['PC1_1']
        del wcs_header['PC2_2']
        try:
            wcs_header['CD2_1'] = wcs_header['PC2_1']
            wcs_header['CD1_2'] = wcs_header['PC1_2']
            del wcs_header['PC2_1']
            del wcs_header['PC1_2']
        except KeyError:
            pass
        h.update(wcs_header)
        return h
        
    def slice(self, cam):
        # Currently using LMI Slice as Placeholder
        # C0_SLICE = np.s_[:, :]
        # C1_SLICE = np.s_[:, :]
        # Crop for data before 2025/11/3
        array_size = 925
        C0_xstart = 48
        C0_ystart = 39
        C1_xstart = 45
        C1_ystart = 36
        C0_SLICE = np.s_[C0_ystart:C0_ystart+array_size, C0_xstart:C0_xstart+array_size]
        C1_SLICE = np.s_[C1_ystart:C1_ystart+array_size, C1_xstart:C1_xstart+array_size]

        C0_SLICE = np.s_[:, :]  # slice is already performed by the detector software
        C1_SLICE = np.s_[:, :]  # slice is already performed by the detector software

        slicedict = {'C0': C0_SLICE, 'C1': C1_SLICE}
        return slicedict[cam]
        
    def get_cam_sat(self, h, idx):
        sat = float(h['SATURATE'])

        return sat
    
    def get_cam_gain(self, h, idx):
        gain = h['GAIN']
        
        return gain      

    def get_exptime(self, h):
        return h['EXPTIME']

    def get_filter(self, h, cam):
        cam_filter = h['FILTER{}'.format(h['ASICINDX']+1)].strip()
        filter2 = h['FILTER3'].lower().strip()
        filter3 = h['FILTER4'].lower().strip()
        if filter2 == 'blank':
            h['FILTER'] = 'blank'
        elif filter2 != 'open':
            if cam_filter.lower().strip() == 'open':
                h['FILTER'] = filter2
            else:
                h['FILTER'] = "{}+{}".format(filter2.upper(), cam_filter)
        else:
            h['FILTER'] = cam_filter
        print(h['FILTER'])
        return h['FILTER'].strip()
            
    def get_centered_filter(self, h, idx):
        return h['FILTER']  

    def change_file_names(self, files):
        obstype_postdict = {
            'SKY_FLAT': 'f', 'DOME_FLAT': 'f',
            'DARK': 'd', 'SCIENCE': 'o', 'SCIENCE CROWDED': 'o'}

        for _file in files:
            pyim = pf.open(_file)
            h = pyim[0].header

            obstype_post = obstype_postdict[h['OBJTYPE'].upper().strip()]
            idate = h['DATE-BEG']
            camnum = h['ASICINDX']
            date_format = '%Y-%m-%dT%H:%M:%S.%f'
            new_date_format = '%Y%m%dT%H%M%S'
            newname = dt.strptime(idate, date_format).strftime(new_date_format) + 'C' + str(camnum)
            newname += obstype_post + '.fits'
            print(newname)
            os.rename(_file, newname)

        return
        
    def original_file_format(self):
        file_format = '????????.rimas.????.??.fits'
        return file_format  


instrument_dict = {'ratir': ratir(), 'lmi': lmi(), 'rimas': rimas()}
