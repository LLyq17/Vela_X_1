# -*- coding: utf-8 -*-
"""
  Author:  jilong 
  Date: 28/02/2020
  Purpose: this is a script to provide the pulse-resolved spectral analysis for Insight-HXMT data
  25/03/2020 fix a bug caused by time_ref
  11/04/2020 fix a bug ceil -> floor 
  12/04/2020 fix a bug, the gti of ME bkg MUST BE RECREATED
  19/04/2020 如果不产生新的bkg，在合并能谱的时候会引发错误的曝光时间
"""
import astropy.io.fits as pyfits
import os, glob
import numpy as np
import argparse
import stingray.pulse as STpulse
import stingray.gti as STgti
from scipy import interpolate 
import binaryCorr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('text', usetex=True)
#rc('font', family='serif')

import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13

#-------------------
#    colored print
HEADER = '\033[95m'
OKBLUE = '\033[94m'
OKGREEN = '\033[92m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'
BOLD = "\033[1m"
#------------------

def create_gti(evt, time_col, binaryPar, freq, phase_bin, time_ref, inst):
    fits = pyfits.open(evt)
    PI = fits[1].data['PI']
    time_obs = fits[1].data['time']  # without the barycentric correction
    time = fits[1].data[time_col]  # with the barycentric correction
    gti = np.array(pyfits.open(evt)[2].data.tolist())
    # ---------------binary correction--------------------
    if binaryPar['pb'] != 0:
        tor = np.empty_like(time)
        for i in range(len(time)):
            tor[i] = binaryCorr.Cor(time[i] / 86400. + hxmt_MJD, binaryPar['t0'], binaryPar['e'], binaryPar['pb'], binaryPar['pbdot'], binaryPar['omega'], binaryPar['omdot'], binaryPar['asini'], binaryPar['gamma'])
        time = time + tor
    # -----------------------------------
    Tobs2Tpulsar = interpolate.interp1d(time_obs, time, fill_value='extrapolate', kind='linear')
    Tpulsar2Tobs = interpolate.interp1d(time, time_obs, fill_value='extrapolate', kind='linear')
    gti = Tobs2Tpulsar(gti)  # GTI in Pulsar Time sys
    # -----------------------------------
    time = time - time_ref
    print(time)
    print(len(time))
    gti = gti - time_ref
    # -----------------------------------
    # calculate how many cycles do we have
    def phase(t, freq):
        #return freq * t + 0.5 * freq_dot * t**2
        return STpulse.pulse_phase(t, *freq, to_1=False)
    #time_steps = np.linspace(time.min(), time.max(), 10000)
    TimeToPhase = interpolate.interp1d(time, phase(time, freq), fill_value='extrapolate', kind='linear')
    PhaseToTime = interpolate.interp1d(phase(time, freq), time, fill_value='extrapolate', kind='linear')
    #print("#####" * 100)
    
    loopMin = np.int64(np.floor(phase(time.min(), freq)))
    loopMax = np.int64(np.ceil(phase(time.max(), freq)))
    Nloop = loopMax - loopMin
    #print(loopMin, loopMax)
    # --------------------------------------------------------------------------
    # get the time in Pulsar Time sys
    def get_time(phase_bin, Nloop):
        p0, p1 = phase_bin
        t0 = np.array([PhaseToTime(p0 + i + loopMin) for i in range(Nloop)])
        t1 = np.array([PhaseToTime(p1 + i + loopMin) for i in range(Nloop)])
        return np.vstack([t0, t1])

    phase_bin = np.array(phase_bin)
    if phase_bin.ndim == 1:
        time_out = get_time(phase_bin, Nloop)
    else:
        time_tmp = []
        for t in range(phase_bin.shape[0]):
            time_tmp.append(get_time(phase_bin[t], Nloop))
        time_out = np.hstack(time_tmp)
        
    mk = np.argsort(time_out[0])
    time_out = np.vstack([time_out[0][mk], time_out[1][mk]])


    # to make it clear I also provide a plot
    if inst == 'LE':
        mk = (PI > 106) & (PI < 1169)
    if inst == "ME":
        mk = (PI > 86) & (PI < 460)
    if inst == "HE":
        mk = (PI > 7) & (PI < 44)
    plt.subplots_adjust(top=0.95, bottom=0.12, left=0.15, right=0.95, hspace=0.01)
    phase_bins, profile, profile_err = STpulse.fold_events(time[mk], *freq, nbin=50, gtis=gti, expocorr=True)
    plt.errorbar(np.hstack([phase_bins, phase_bins + 1]), np.hstack([profile, profile]), xerr=1 / 100., yerr=np.hstack([profile_err, profile_err]), drawstyle='steps-mid')
    plt.xlim(0, 2)
    plt.xlabel("Phase", fontsize=13)
    #plt.show()
    if inst == 'LE':
        plt.ylabel(r"LE (1-10 keV)", fontsize=13)
    if inst == "ME":
        plt.ylabel(r"ME (8-30 keV)$", fontsize=13)
    if inst == "HE":
        plt.ylabel(r"HE (25-100 keV)", fontsize=13)

    if phase_bin.ndim == 1:
        plt.fill_betweenx([(profile - profile_err).min(), (profile.max() + profile_err).max()], phase_bin[0], phase_bin[1], facecolor='red', alpha=0.2)
    else:
        for t in range(phase_bin.shape[0]):
            plt.fill_betweenx([(profile - profile_err).min(), (profile.max() + profile_err).max()], phase_bin[t][0], phase_bin[t][1], facecolor='red', alpha=0.2)

    plt.ylim([(profile - profile_err).min(), (profile.max() + profile_err).max()])
    plt.savefig("save_tmp.pdf")
    plt.cla()

    #print(gti.shape, time_out.shape)
    # cross_two_gti = STgti.cross_two_gtis(gti, time_out.T) #it is very slow if the freq is high and almost useless
    #print("!!!")
    # exposure = (cross_two_gti[:, 1] - cross_two_gti[:, 0]).sum()  # the exposure in different time sys is very slightly different, here I close this step
    # return Tpulsar2Tobs(cross_two_gti + time_ref), exposure
    
    
    # ---------------------------------------------
    # ---------------------------------------------
    
    return Tpulsar2Tobs(time_out.T + time_ref), 0

# -------------------------------------------------------------------------------
def lebkgmap(blind_det, gtifile, lc_spec_name, chmin, chmax, outnam_prefix):
    TPYE = 'spec'
    string = "lebkgmap %s %s %s %s %s %s %s" % (TPYE, blind_det, gtifile, lc_spec_name, chmin, chmax, outnam_prefix)
    os.system(string)

def lespecgen(lescreenEVT, pathout, outfile):
    string = "lespecgen evtfile=%s outfile=%s/%s_LE userdetid='%s'  eventtype=1  clobber=yes" % (lescreenEVT, pathout, outfile, ledetidsmall)
    os.system(string)

def lerspgen(phafile, outfile, attfile, tempfile):
    string = "lerspgen  phafile=%s outfile=%s attfile=%s tempfile=%s  clobber=yes  ra=-1  dec=-91" % (phafile, outfile, attfile, tempfile)
    os.system(string)
# -------------------------------------------------------------------------------
def hxmtscreen(evt, out, gti):
    com = 'hxmtscreen evtfile=%s outfile=%s usergtifile=%s userdetid="NONE"' % (evt, out, gti)
    os.system(com)    

# -------------------------------------------------------------------------------
def write_gti_fits(gti):
    col1 = pyfits.Column(name='start', format='D', array=gti[:, 0])
    col2 = pyfits.Column(name='stop', format='D', array=gti[:, 1])
    cols = pyfits.ColDefs([col1, col2])
    hdu = pyfits.BinTableHDU.from_columns(cols)
    hdu.writeto('gti_tmp.fits', overwrite=True)
    
# -------------------------------------------------------------------------------
def grppha(infil, outfil, expr):
    string = "grppha '%s'  '!%s' '%s'  exit" % (infil, outfil, expr)
    os.system(string)

# -------------------------------------------------------------------------------
def mespecgen(mescreenEVT, pathout, outfile, medead):
    string = "mespecgen evtfile=%s outfile=%s/%s_ME deadfile=%s  userdetid='%s' clobber=yes" % (mescreenEVT, pathout, outfile, medead, medetidsmall)
    os.system(string)

def mebkgmap(blind_det, ehkfile, gtifile, deadtime, temperature, lc_spec_name, chmin, chmax, outnam_prefix, det_status):
    TPYE = "spec"
    string = "mebkgmap %s %s %s %s %s %s %s %s %s  %s  %s" % (TPYE, blind_det, ehkfile, gtifile, deadtime, temperature, lc_spec_name, chmin, chmax, outnam_prefix, det_status)
    os.system(string)


def merspgen(phafile, outfile, attfile):
    string = "merspgen  phafile=%s outfile=%s attfile=%s  clobber=yes ra=-1 dec=-91" % (phafile, outfile, attfile)
    os.system(string)

# -------------------------------------------------------------------------------
def hespecgen(hescreenEVT, pathout, outfile, hedead, userdetid):
    string = 'hespecgen evtfile=%s outfile=%s/%s_HE deadfile=%s  userdetid="%s"  clobber=yes' % (hescreenEVT, pathout, outfile, hedead, userdetid)
    os.system(string)
    
def hebkgmap(blind_det, ehkfile, gtifile, deadtime, lc_spec_name, chmin, chmax, outnam_prefix):
    TPYE = 'spec'
    string = "hebkgmap %s %s %s %s %s %s %s %s %s" % (TPYE, blind_det, ehkfile, gtifile, deadtime, lc_spec_name, chmin, chmax, outnam_prefix)
    os.system(string)

def herspgen(phafile, outfile, attfile):
    string = "herspgen phafile=%s outfile=%s attfile=%s ra=-1 dec=-91  clobber=yes" % (phafile, outfile, attfile)
    os.system(string)
def hhe_spec2pi(src_lst, bkg_lst, rsp_lst, src_pi, bkg_pi, rsp_rsp):
    string = "hhe_spec2pi  %s   %s   %s   %s   %s  %s" % (src_lst, bkg_lst, rsp_lst, src_pi, bkg_pi, rsp_rsp)
    os.system(string)
    
def writetofit_MEnewBad(new_bad_det_list, new_bad):
    if os.path.exists(new_bad):
        os.system("rm %s" % new_bad)
    hdu = pyfits.PrimaryHDU()
    c1 = pyfits.Column(name='DetID', array=new_bad_det_list, format='I')
    c2 = pyfits.Column(name='TIMERANGE', array=[0 for i in range(len(new_bad_det_list))], format='20A')
    c3 = pyfits.Column(name='TIMERANGE2', array=['INDEF' for i in range(len(new_bad_det_list))], format='20A')
    c4 = pyfits.Column(name='TYPE', array=['bad' for i in range(len(new_bad_det_list))], format='20A')
    c5 = pyfits.Column(name='STATUS', array=[0 for i in range(len(new_bad_det_list))], format='B')
    hdu2 = pyfits.BinTableHDU.from_columns([c1, c2, c3, c4, c5])
    hdu2.header['EXTNAME'] = 'detectorStatus'
    hdu2.header['TELESCOP'] = "HXMT"
    hdu2.header['INSTRUME'] = "ME"
    hdu2.header['VERSION'] = "1.00"
    hdu2.header['CREATOR'] = "JL"
    hdu2.header['UNITNUM'] = "1728"
    hdul = pyfits.HDUList([hdu, hdu2])
    hdul.writeto(new_bad)
    
hxmt_MJD = 55927.00076601852
if __name__ == "__main__":
    readme = '''
    This is a script to perform a phase-resolved analysis in a convenient way for Insight-HXMT users. Before running this script, screened event files that only include good time intervals you defined, as well as the corresponding auxiliary files, should be prepared. It is recommended to get these files by using the pipeline "HXMTanalysis.py" to reduce some pains. However, this script is also compatible with other pipelines when you use the "-NotFromHXMTpipine" flag and specify relevant files by using options like "-LEscreen". In order to calculate the phase, a trial "frequency array ([F0, F1, F2, ...])" should be given. You can also define the reference time (-time_ref) for phase=0. By default, the reference time is the starting of LE or ME or HE that will be shown after running. 
    If your source is in a binary system, the binary movement should be considered. In this case, the ephemeris is needed. Here is an example to extract a spectrum from the phase 0.1-0.2&&0.7-0.8, assuming that the frequency is 1Hz, the frequency_dot is 0.001, and the prefix name of outfiles is phase_0.1-0.2and0.7-0.8. 
    ---------------
    python phase_resolved_spectroscopy.py outpath phase_0.1-0.2and0.7-0.8  -freq '1, 0.001'  -pipineOutput P021400100101-20190120-01-01 -phase_bin_min 0.1 0.7 -phase_bin_max 0.2 0.8 -binPar_pb 9599040.0 -binPar_t0  53613 -binPar_asini 267 -binPar_e 0.47 -binPar_omega 2.27
    ---------------
    The screened events are stored in the folder "P021400100101-20190120-01-01" which was produced by the "HXMTanalysis.py" script previously. The -phase_bin_min and -phase_bin_max arguments are used to define the phase intervals. The arguments starting with "-binPar_" are the binary parameters.
    
    After running this script, the phase-resolved spectrum will be found in "outpath", accompanied with a background file and a response file and a PDF file which contains the pulse profile of this source and the phase bins we are interested in are highlighted. 

    '''
    
    parser = argparse.ArgumentParser(description=readme, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("outpath", help="the path of outfiles", type=str)
    parser.add_argument("outfile_name", help="the prefix of outfiles", type=str)
    parser.add_argument("-freq", help="the frequency array, such as \"-freq 'F0, F1, F2,...'\"", type=str)
    #parser.add_argument("-fdot", help="the frequency derivative at the reference time", type=float, default=0)
    # ------------------------------------
    parser.add_argument("-phase_bin_min", nargs='+', type=float, default=[0], required=True,
                        help="a list of starting points of phase bins that much be between 0 and 1. Note that NO overlap is allowed!")
    parser.add_argument("-phase_bin_max", nargs='+', type=float, default=[0.9999], required=True,
                        help="a list of ending points of phase bins that much be between 0 and 1")
    parser.add_argument("-NotFromHXMTpipine", action="store_true", default=False,
                        help="use this flag only if you have to specify the event files you want to use and the corresponding auxiliary files")
    parser.add_argument("-pipineOutput", type=str,
                        help="where you save the outputs when runing HXMTanalysis.py, this is a recommended way if you have done your data analysis with the script HXMTanalysis.py")
    # ------------------------------------

    parser.add_argument("-time_ref", type=float, default=0,
                        help="(optional) unit: s, the reference time for the timing analysis (in Pulsar sys assuming HMXT's MJD); a nagetive value means the Tstart of this obsID")
    # -----------------------------------
    parser.add_argument("-binPar_pb", help="binary orbital period (s)", type=float, default=0)
    parser.add_argument("-binPar_asini", help="projected semi-major axis (light seconds)", type=float, default=0)
    parser.add_argument("-binPar_e", help="orbital eccentricity", type=float, default=0)
    parser.add_argument("-binPar_t0", help="barycentric time (in MJD(TDB)) of periastron ", type=float, default=0)
    parser.add_argument("-binPar_omega", help="longitude of periastron (radian degrees)", type=float, default=0)
    parser.add_argument("-binPar_omdot", help="first derivative of omdot (degrees per Julian year)", type=float, default=0)
    parser.add_argument("-binPar_gamma", help="time-dilation and gravitational redshift parameter (s)", type=float, default=0)
    parser.add_argument("-binPar_pbdot", help="first derivative of pb", type=float, default=0)

    #parser.add_argument("-bkgFlag", action="store_true", default=False,
                        #help="force to run the phase-resolved bkg estimation")

    parser.add_argument("-time_col", type=str, default='TDB',
                        help="the time column to be used in the timing analysis")
    # -----------------------------------
    parser.add_argument("-LEscreen", help="a LE event file used to extract the spectrum", type=str)
    parser.add_argument("-MEscreen", help="a ME event file used to extract the spectrum", type=str)
    parser.add_argument("-HEscreen", help="a HE event file used to extract the spectrum", type=str)

    parser.add_argument("-LEblind", help="a LE event file that contains blind detectors, maybe the same as -LEscreen", type=str)
    parser.add_argument("-MEblind", help="a ME event file that contains blind detectors, maybe the same as -MEscreen", type=str)
    parser.add_argument("-HEblind", help="a HE event file that contains blind detectors, maybe the same as -HEscreen", type=str)

    parser.add_argument("-ehk", help="the ehk file", type=str)
    parser.add_argument("-leth", help="the LE temperature file", type=str)
    parser.add_argument("-att", help="the attitude file", type=str)
    parser.add_argument("-medead", help="the ME deadtime file", type=str)
    parser.add_argument("-meth", help="the LE temperature file", type=str)
    parser.add_argument("-meDetectorStatus", help="the status of ME pixels which is created by using megti", type=str, default="")
    parser.add_argument("-hedead", help="the HE deadtime file", type=str)
    
    parser.add_argument("-LeDetID", help="an ascii file including user-defined detectors of LE", default='', type=str)
    parser.add_argument("-MeDetID", help="an ascii file including user-defined detectors of ME", default='', type=str)
    parser.add_argument("-HeDetID", help="an ascii file including user-defined detectors of HE", default='', type=str)
    

    args = parser.parse_args()

    if args.NotFromHXMTpipine == False:
        le_evt = '%s/LE/lescreen_v204.fits' % args.pipineOutput
        me_evt = '%s/ME/mescreen_v204.fits' % args.pipineOutput
        he_evt = '%s/HE/hescreen_v204.fits' % args.pipineOutput

        #le_blind = '%s/LE/lescreen_v204.fits.blind' % args.pipineOutput
        #me_blind = '%s/ME/mescreen_v204.fits.blind' % args.pipineOutput
        #he_blind = '%s/HE/hescreen_v204.fits.blind' % args.pipineOutput
       
        le_blind = '%s/LE/lescreen_v204.fits' % args.pipineOutput
        me_blind = '%s/ME/mescreen_v204.fits' % args.pipineOutput
        he_blind = '%s/HE/hescreen_v204.fits' % args.pipineOutput

        # ----------------------------------
        ehkfile = '%s/ehk.fits' % args.pipineOutput
        leth = "%s/LE/leth.fits" % args.pipineOutput
        att = "%s/att.fits" % args.pipineOutput
        medead = '%s/ME/deadfile.fits' % args.pipineOutput
        meth = '%s/ME/meth.fits' % args.pipineOutput
        me_det_status = '%s/ME/newmedetectorstatus.fits' % args.pipineOutput
        hedead = '%s/HE/deadtime.fits' % args.pipineOutput
        # -----------------------------------
    else:
        le_evt = args.LEscreen
        me_evt = args.MEscreen
        he_evt = args.HEscreen

        le_blind = args.LEblind
        me_blind = args.MEblind
        he_blind = args.HEblind
        # -----------------------
        ehkfile = args.ehk
        leth = args.leth
        att = args.att
        medead = args.medead
        meth = args.meth
        me_det_status = args.meDetectorStatus
        hedead = args.hedead

    outpath = args.outpath
    outfile = args.outfile_name

    freq = [float(item) for item in args.freq.split(',')]
    time_col = args.time_col

    phase_bin = np.vstack([args.phase_bin_min, args.phase_bin_max]).T
    binaryPar = {
        'pb': args.binPar_pb,  # /* Orbital period (s) */
        'asini': args.binPar_asini,  # /* Projected semi-major axis (light seconds) */
        'e': args.binPar_e,  # /* orbital eccentricity */
        't0': args.binPar_t0,  # /* Barycentric time (in MJD(TDB)) of periastron */
        'omega': args.binPar_omega,  # /* Longitude of periastron (radian degrees) */
        'omdot': args.binPar_omdot,  # /* First derivative of omdot (degrees per Julian year) */
        'gamma': args.binPar_gamma,  # /* Time-dilation and gravitational redshift parameter (s) */
        'pbdot': args.binPar_pbdot,  # /* First derivative of pb */
    }

    #if freq > 0.1:
        #bkgFlag = False
    #else:
        #bkgFlag = True
        
    #if args.bkgFlag == True:
        #bkgFlag = True

    # -------------------

    if args.LeDetID != '':
        ledetidsmall = np.genfromtxt(args.LeDetID, str, delimiter='|').reshape(-1)[0]
    else:
        ledetidsmall = "0, 2, 3, 4, 6, 7, 8, 9, 10, 12, 14, 20, 22, 23, 24, 25, 26, 28, 30, 32, 34, 35, 36, 38, 39, 40, 41, 42, 44, 46, 52, 54, 55, 56, 57, 58, 60, 61, 62, 64, 66, 67, 68, 70, 71, 72, 73, 74, 76, 78, 84, 86, 88, 89, 90, 92, 93, 94"
    #  modified by qi---remove  29, 87 ID of LE 
    if args.MeDetID != '':
        medetidsmall = np.genfromtxt(args.MeDetID, str, delimiter='|').reshape(-1)[0]
    else:
        medetidsmall = "0-7 11-17 18-25 29-35 36-43 47-53"
        

    if args.HeDetID != '':
        hedetid = np.genfromtxt(args.HeDetID, str, delimiter='|').reshape(-1)[0]
    else:
        hedetid = "0-15,17"

    # -------------------    
    #print("###" * 10, freq)
    # -----------start to do analysis----------------
    if args.time_ref >= 0:
        time_ref = args.time_ref 
    else:
        time_ref_tmp = []
        for evt in [le_evt, me_evt, he_evt]:
            if os.path.exists(evt):
                time_ref_tmp.append(pyfits.open(evt)[1].header['tstart'])
            else:
                print("%s does not exist" % evt)
        time_ref = np.min(time_ref_tmp)
    
    # ---------------------------
    if not os.path.exists(outpath):
        os.system("mkdir -p %s" % outpath)
    # ---------------------------
    # for evt in [le_evt, me_evt, he_evt]:
    for evt in [le_evt, me_evt, he_evt]:
        if os.path.exists(evt):
            inst = pyfits.open(evt)[1].header['INSTRUME'].strip().upper()
            if inst == "LE":
                gti_tmp, exposure_tmp = create_gti(evt, time_col, binaryPar, freq, phase_bin, time_ref, inst)
                #print(gti_tmp)
                # ---
                gti_old = np.array(pyfits.open(evt)[2].data.tolist())
                gti_tmp = STgti.cross_two_gtis(gti_old, gti_tmp)
                # ---
                np.savetxt("gti_tmp.txt", gti_tmp)
                # --------------plot the pulse profile-----------------------
                os.system("mv save_tmp.pdf  %s/%s_LE.pdf" % (outpath, outfile))
                # -----------------------------------------
                hxmtscreen(evt, "evt_tmp.fits", 'gti_tmp.txt')
                # -----------------------------------------
                lespecgen("evt_tmp.fits", outpath, outfile)

                glob_list = glob.glob('%s/%s_LE_g*.pha' % (outpath, outfile))
                LE_spec_file = sorted(glob_list, key=os.path.getctime)[-1]
                os.system("mv   %s   %s/%s_LE.pha" % (LE_spec_file, outpath, outfile))
                
                # ------------------------------------------
                # have to make a new LE gti in fits format ??
                write_gti_fits(gti_tmp)
                # -------------------------------------------
                f = open("lebkgmap.tmp", "w")
                f.write("%s/%s_LE.pha\n" % (outpath, outfile))
                f.close()
                if True:
                    #hxmtscreen(le_blind, "evt_blind_tmp.fits", 'gti_tmp.txt')
                    lebkgmap("evt_tmp.fits", "gti_tmp.fits", 'lebkgmap.tmp', 0, 1535, "%s/%s_LE_bkg" % (outpath, outfile))  # here I use le_blind, check with Mingyu
                    # modified by qi  106, 1169 --> 0, 1535
                    #os.system("rm evt_blind_tmp.fits")
                    # --rescale---
                    if os.path.exists("%s/%s_LE_bkg.pha" % (outpath, outfile)):
                        
                        #glob_list = glob.glob('%s/lespec_g0_0-94.pha' % args.pipineOutput) #modified by qi
                        #LE_spec_old = sorted(glob_list, key=os.path.getctime)[0] #removed by qi 
                        LE_spec_old = '%s/lespec_g0_0-94.pha' % args.pipineOutput 

                        print('####' * 10, LE_spec_old)
                        det_nu_1 = len(pyfits.open(LE_spec_old)[-1].data['ID']) * 1.0
                        det_nu_2 = len(pyfits.open("%s/%s_LE.pha" % (outpath, outfile))[-1].data['ID']) * 1.0
                        # ---
                        fits = pyfits.open("%s/%s_LE_bkg.pha" % (outpath, outfile), 'update')
                        if "rescale" not in list(fits[1].header):
                            fits[1].header['exposure'] = fits[1].header['exposure'] * det_nu_1 / det_nu_2
                            fits[1].header['rescale'] = det_nu_2 / det_nu_1
                            fits.flush()
                        fits.close()                    
                # -------------------------------------------
                lerspgen("%s/%s_LE.pha" % (outpath, outfile), "%s/%s_LE.rsp" % (outpath, outfile), att, leth)
                
                grppha("%s/%s_LE.pha" % (outpath, outfile), "%s/%s_LE.pha" % (outpath, outfile), "group min 30")
                grppha("%s/%s_LE.pha" % (outpath, outfile), "%s/%s_LE.pha" % (outpath, outfile), "chkey respfile  %s_LE.rsp" % outfile)
                
                if True:
                    grppha("%s/%s_LE.pha" % (outpath, outfile), "%s/%s_LE.pha" % (outpath, outfile), "chkey backfile  %s_LE_bkg.pha" % outfile)
                
                # --------Cleaning----------------------------------      
                os.system("rm gti_tmp.fits  lebkgmap.tmp  gti_tmp.txt  evt_tmp.fits")
            elif inst == "ME":
                gti_tmp, exposure_tmp = create_gti(evt, time_col, binaryPar, freq, phase_bin, time_ref, inst)
                # ---
                gti_old = np.array(pyfits.open(evt)[2].data.tolist())
                gti_tmp = STgti.cross_two_gtis(gti_old, gti_tmp)
                # ---                
                np.savetxt("gti_tmp.txt", gti_tmp)
                write_gti_fits(gti_tmp)
                # --------------plot the pulse profile-----------------------
                os.system("mv save_tmp.pdf  %s/%s_ME.pdf" % (outpath, outfile))
                # -----------------------------------------
                hxmtscreen(evt, "evt_tmp.fits", 'gti_tmp.txt')
                # -----------------------------------------
                mespecgen("evt_tmp.fits", outpath, outfile, medead)

                glob_list = glob.glob('%s/%s_ME_g*.pha' % (outpath, outfile))
                ME_spec_file = sorted(glob_list, key=os.path.getctime)[-1]

                os.system("mv   %s   %s/%s_ME.pha" % (ME_spec_file, outpath, outfile))
                # -----------------------------------------
                write_gti_fits(gti_tmp)
                f = open("mebkgmap.tmp", "w")
                f.write("%s/%s_ME.pha\n" % (outpath, outfile))
                f.close()
                if True:
                    #glob_list = glob.glob('%s/mespec_g0_0-53.pha' % args.pipineOutput)  #modified by qi
                    #old_spectr = sorted(glob_list, key=os.path.getctime)[0]
                    old_spectr ='%s/mespec_g0_0-53.pha' % args.pipineOutput

                    spectr = '%s/%s_ME.pha' % (outpath, outfile)
                    old_bad = me_det_status
                    new_bad = './user_medetectorstatus.fits'

                    old_ID = pyfits.open(old_spectr)[-1].data['id']
                    new_ID = pyfits.open(spectr)[-1].data['id']
                
                    mk = np.array([item not in new_ID for item in old_ID])
                    ID = old_ID[mk]
                
                    bad_det_list = []
                    for i in range(len(ID)):
                        bad_det_list.append(ID[i] * 32 + np.arange(0, 32))
                    if len(bad_det_list) > 0:
                        bad_det_list = np.hstack(bad_det_list)
                    else:
                        bad_det_list = np.array([])
                    # ----
                    old_bad_det_list = pyfits.open(old_bad)[1].data['DetID']
                    new_bad_det_list = np.array(list(set(np.hstack([bad_det_list, old_bad_det_list]))))
                    # ---
                    writetofit_MEnewBad(new_bad_det_list, new_bad)
                    # ---                    
                    
                    #hxmtscreen(me_blind, "evt_blind_tmp.fits", 'gti_tmp.txt')                    
                    mebkgmap("evt_tmp.fits", ehkfile, 'gti_tmp.fits', medead, meth, 'mebkgmap.tmp', 0, 1023, "%s/%s_ME_bkg" % (outpath, outfile), new_bad)
                    #os.system("rm evt_blind_tmp.fits")                   
                # -----------------------------------------
                merspgen("%s/%s_ME.pha" % (outpath, outfile), "%s/%s_ME.rsp" % (outpath, outfile), att)
                # -----------------------------------------                                
                grppha("%s/%s_ME.pha" % (outpath, outfile), "%s/%s_ME.pha" % (outpath, outfile), "group min 30")
                if True:
                    grppha("%s/%s_ME.pha" % (outpath, outfile), "%s/%s_ME.pha" % (outpath, outfile), "chkey backfile  %s_ME_bkg.pha" % outfile)
                grppha("%s/%s_ME.pha" % (outpath, outfile), "%s/%s_ME.pha" % (outpath, outfile), "chkey respfile  %s_ME.rsp" % outfile)
                # -----------------------------------------
                os.system("rm gti_tmp.fits mebkgmap.tmp  gti_tmp.txt  evt_tmp.fits  %s" % new_bad)

            elif inst == "HE":
                gti_tmp, exposure_tmp = create_gti(evt, time_col, binaryPar, freq, phase_bin, time_ref, inst)
                # ---
                gti_old = np.array(pyfits.open(evt)[2].data.tolist())
                gti_tmp = STgti.cross_two_gtis(gti_old, gti_tmp)
                # ---                
                np.savetxt("gti_tmp.txt", gti_tmp)
                write_gti_fits(gti_tmp)
                # --------------plot the pulse profile-----------------------
                os.system("mv save_tmp.pdf  %s/%s_HE.pdf" % (outpath, outfile))
                # -----------------------------------------
                hxmtscreen(evt, "evt_tmp.fits", 'gti_tmp.txt')
                # -----------------------------------------   
                #for i in range(18):
                #    hespecgen("evt_tmp.fits", outpath, outfile, hedead, i)
                #    if i != 16:
                #        herspgen("%s/%s_HE_g0_%d.pha" % (outpath, outfile, i), "%s/%s_HE_%d.rsp" % (outpath, outfile, i), att)
                hespecgen("evt_tmp.fits", outpath, outfile, hedead, hedetid)

                glob_list = glob.glob('%s/%s_HE_g*.pha' % (outpath, outfile))
                HE_spec_file = sorted(glob_list, key=os.path.getctime)[-1]

                os.system("mv   %s   %s/%s_HE.pha" % (HE_spec_file, outpath, outfile))
                herspgen("%s/%s_HE.pha" % (outpath, outfile), "%s/%s_HE.rsp" % (outpath, outfile), att)
                #HE_expo = pyfits.open("%s/%s_HE_g0_0.pha" % (outpath, outfile))[1].header['exposure']
                # ----start to estimate the background-----
                #f = open('he_pha_bkg_tmp.txt', 'w')
                #for i in range(18):
                #    s = "%s/%s_HE_g0_%d.pha" % (outpath, outfile, i)
                #    f.write("%s\n" % s)
                #f.close()
                f = open('he_pha_bkg_tmp.txt', 'w')
                s = "%s/%s_HE.pha" % (outpath, outfile)
                f.write("%s\n" % s)
                f.close()
                # -----------------------------------------
                if True:
                    #hxmtscreen(he_blind, "evt_blind_tmp.fits", 'gti_tmp.txt')
                    hebkgmap("evt_tmp.fits", ehkfile, 'gti_tmp.fits', hedead, 'he_pha_bkg_tmp.txt', 0, 255, "%s/%s_HE_bkg" % (outpath, outfile)) # modified by qi 0, 256 --> 0, 255
                    #os.system("rm evt_blind_tmp.fits")
                else:
                    #-------using the phase averaged spectrum, but the exposure should be rescaled
                    #for i in range(18):
                    #    os.system("cp  %s/HE/*_HE_bkg_%d.pha  %s/%s_HE_bkg_%d.pha" % (args.pipineOutput, i, outpath, outfile, i))
                    #    f = pyfits.open("%s/%s_HE_bkg_%d.pha" % (outpath, outfile, i), 'update')
                    #    expo = f[1].header['exposure']
                    #    f[1].header['exposure'] = (gti_tmp[:, 1] - gti_tmp[:, 0]).sum()  # expo * (phase_bin[:, 1] - phase_bin[:, 0]).sum()
                    #    f.flush()
                    #    f.close()
                    # os.system("cp  %s/HE/*_HE_bkg.pha  %s/%s_HE_bkg.pha" % (args.pipineOutput, outpath, outfile))
                    pass

                grppha("%s/%s_HE.pha" % (outpath, outfile), "%s/%s_HE.pha" % (outpath, outfile), "group min 30")
                grppha("%s/%s_HE.pha" % (outpath, outfile), "%s/%s_HE.pha" % (outpath, outfile), "chkey respfile  %s_HE.rsp" % outfile)
                
                if True:
                    grppha("%s/%s_HE.pha" % (outpath, outfile), "%s/%s_HE.pha" % (outpath, outfile), "chkey backfile  %s_HE_bkg.pha" % outfile)
                #else:
                #    os.system("rm  %s/%s_HE_bkg.pha" % (outpath, outfile))

                # -----------------------------------------    
                os.system("rm  he_pha_bkg_tmp.txt evt_tmp.fits gti_tmp.fits gti_tmp.txt")

        else:
            print("%s does not exist, so no corresponding spectrum will be created" % evt)
    
    print("******" * 3, "Finish", "******" * 3)
    print("The reference time for the timing analysis is %.10lf (in MJD)" % (time_ref / 86400. + hxmt_MJD))



