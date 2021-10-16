# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 11:40:12 2021

@author: xiao
"""

import os
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import xspec as xs
from xspec import *
mpl.rcParams['xtick.labelsize'] = 13

AllModels.initpackage("mymodel","lmodel.dat","/home/lq/heasoft-6.28/Xspec/src/Qifunc")
AllModels.lmod("mymodel","/home/lq/heasoft-6.28/Xspec/src/Qifunc")

path = "/home/lq/hxmt/Vela_X-1/P0201012_v204/P0201012149/P020101214906-20200323-02-01"
path_out = path + '/'

model_par = ["constant*phabs*(pow*fdcut*cyclabs*cyclabs + gaussian),", "nH,", "E2cyc,", "Width2,", "Depth2,", "E1cyc,", "Width1,", "Depth1,", "PhoIndex,", "pl norm,", "Ecut,", "Efold,", "E_Fe,", "gau sigma,", "gau norm,", "chi square dof,", "const1_1,", "const1_2,", "const1_3"]

for root, dirs, files in os.walk(path):
    print(root,dirs,files)
    if (("hespec_g0_0-17.pha" in files) and ("mespec_g0_0-53.pha" in files) and ("lespec_g0_0-94.pha" in files)):
        print(root)
        os.chdir(root)
        
        fig = plt.figure(figsize=(12,9))
        gs = gridspec.GridSpec(5,1,hspace=0,wspace=0,height_ratios=[3,1,1,1,1])
        fig.subplots_adjust(wspace=0, hspace=0)
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[1,0])
        ax3 = plt.subplot(gs[2,0])
        ax4 = plt.subplot(gs[3,0])
        ax5 = plt.subplot(gs[4,0])
        
        xs.Plot.device='/null'
        xs.Plot.xAxis="keV"
        
        xs.AllData("1:1 lespec_g0_0-94.pha 2:2 mespec_g0_0-53.pha 3:3 hespec_g0_0-17.pha")
        
        xs.AllData(1).ignore("**-2. 10.-**")
        xs.AllData(2).ignore("**-9. 29.-**")
        xs.AllData(3).ignore("**-27. 105.-**")
        
        model_string = "constant*phabs*(pow*fdcut*cyclabs*cyclabs + gaussian)"
        m1 = xs.Model("%s"%model_string)
        m2 = xs.AllModels(2)
        m3 = xs.AllModels(3)

        m1.setPars({8:"45., .01, 35, 35, 55, 55"},{9:9.3},{13:"25., .01, 18, 18, 30, 30"},{14:3.9})
        m2(1).link = ""
        m3(1).link = ""
        m1(1).frozen = True
        m2.setPars("1., .001, .85, .85, 1.05, 1.05")
        m3.setPars("1., .001, .85, .85, 1.05, 1.05")
        xs.AllModels.show()

        #if '102' in root:
            #m1.setPars("1.,-1","1.4,-1",'90.0,,40.,40.,120.,120.','15.,,14,14,16,16', 50.,0.4,{8:12.},{9:17.},{10:2.},{12:6.5},{13:0.1})
        #elif '313' in root:
            #m1.setPars("1.,-1","1.4,-1",'90.0,,40.,40.,120.,120.','11.,,10,10,12,12', 50.,0.4,{8:12.},{9:17.},{10:2.},{12:6.5},{13:0.1})
        #else:
            #m1.setPars("1.,-1","1.4,-1",'90.0,,40.,40.,120.,120.',12., 50.,0.4,{8:12.},{9:17.},{10:2.},{12:6.5},{13:0.1})
        
        xs.Fit.renorm()
        xs.Fit.nIterations = 100000
        xs.Fit.query = "yes"
        xs.Fit.perform()
        print("model=%s, chi^2 = "%model_string,xs.Fit.statistic,"/",xs.Fit.dof,"=",xs.Fit.statistic/xs.Fit.dof)
        
        temp_const_m1 = m1.constant.factor.values[0]
        temp_phabs_nH = m1.phabs.nH.values[0]
        temp_cyclabs_E2 = m1(8).values[0]
        temp_cyclabs_W2 = m1(9).values[0]
        temp_cyclabs_D2 = m1(7).values[0]
        temp_cyclabs_E1 = m1(13).values[0]
        temp_cyclabs_W1 = m1(14).values[0]
        temp_cyclabs_D1 = m1(12).values[0]
        temp_pl_gamma = m1.powerlaw.PhoIndex.values[0]
        temp_pl_norm = m1.powerlaw.norm.values[0]
        temp_fdcut_cutoffE = m1.fdcut.cutoffE.values[0]
        temp_fdcut_foldE = m1.fdcut.foldE.values[0]
        temp_gaussian_LineE = m1.gaussian.LineE.values[0]
        temp_gaussian_Sigma = m1.gaussian.Sigma.values[0]
        temp_gaussian_norm = m1.gaussian.norm.values[0]
        temp_const_m2 = m2.constant.factor.values[0]
        temp_const_m3 = m3.constant.factor.values[0]
        
        model_par = np.row_stack((model_par, ["%s%s,"%(root[-3:-1],root[-1]), "%s,"%temp_phabs_nH, "%s,"%temp_cyclabs_E2, "%s,"%temp_cyclabs_W2, "%s,"%temp_cyclabs_D2, "%s,"%temp_cyclabs_E1, "%s,"%temp_cyclabs_W1, "%s,"%temp_cyclabs_D1, "%s,"%temp_pl_gamma, "%s,"%temp_pl_norm, "%s,"%temp_fdcut_cutoffE, "%s,"%temp_fdcut_foldE, "%s,"%temp_gaussian_LineE, "%s,"%temp_gaussian_Sigma, "%s,"%temp_gaussian_norm, "%s,"%(xs.Fit.statistic/xs.Fit.dof), "%s,"%temp_const_m1, "%s,"%temp_const_m2, "%s"%temp_const_m3]))
        
        xs.Plot.setRebin(18, 18, 1)
        xs.Plot.setRebin(29, 29, 2)
        xs.Plot.setRebin(1, 1, 3)
        xs.Plot("uf")
        LEenergy = xs.Plot.x(1)
        LExerror = xs.Plot.xErr(1)
        LEkeVsq = xs.Plot.y(1)
        LEyerror = xs.Plot.yErr(1)
        LEfolded = xs.Plot.model(1)
        
        MEenergy = xs.Plot.x(2)
        MEkeVsq = xs.Plot.y(2)
        MExerror = xs.Plot.xErr(2)
        MEyerror = xs.Plot.yErr(2)
        MEfolded = xs.Plot.model(2)
        
        HEenergy = xs.Plot.x(3)
        HEkeVsq = xs.Plot.y(3)
        HExerror = xs.Plot.xErr(3)
        HEyerror = xs.Plot.yErr(3)
        HEfolded = xs.Plot.model(3)
        
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_ylabel('Photons c$\mathregular{m^{-2}}$ $\mathregular{s^{-1}}$ ke$\mathregular{V^{-1}}$)', size=15)
        ax1.set_xlim(2,105)
        
        ax1.errorbar(LEenergy, LEkeVsq, yerr=LEyerror, xerr=LExerror, fmt='+', markersize=0.2, elinewidth=0.65)#,label='%s'%exeObject)
        ax1.plot(LEenergy, LEfolded, 'r', linestyle='-',linewidth=1)

        ax1.errorbar(MEenergy, MEkeVsq, yerr=MEyerror, xerr=MExerror, fmt='+', markersize=0.2, elinewidth=0.65)
        ax1.plot(MEenergy, MEfolded, 'r', linestyle='-',linewidth=1)

        ax1.errorbar(HEenergy, HEkeVsq, yerr=HEyerror, xerr=HExerror, fmt='+', markersize=0.2, elinewidth=0.65)
        ax1.plot(HEenergy, HEfolded, 'r', linestyle='-',linewidth=1)
        #ax1.legend(fontsize=15)
        #ax1.set_xticks([])
        ax1.tick_params(axis='x', direction='in', which='both', top=True, labelbottom=False)
        ax1.tick_params(axis='y', labelsize=13, direction='out', which='both', right=True)
        
        #                ax1.set_yticks([2,5,10,20],[2,5,10,20])
        
        xs.Plot.add = True
        xs.Plot("model")
        LEenergy_all = xs.Plot.x(1)
        for i in range(len(LEenergy_all)):
            if LEenergy_all[i] >= 10:
                LEhigh_step = i
                break
        #LEmodel = xs.Plot.model(1)
        LE1_all = xs.Plot.addComp(1,1)
        LE2_all = xs.Plot.addComp(2,1)
        LEenergy = LEenergy_all[:LEhigh_step]
        LE1 = LE1_all[:LEhigh_step]
        LE2 = LE2_all[:LEhigh_step]
        
        MEenergy_all = xs.Plot.x(2)
        for i in range(len(MEenergy_all)):
            if MEenergy_all[i] >= 9:
                MElow_step = i
                break
        for i in range(len(MEenergy_all)):
            if MEenergy_all[i] >= 29:
                MEhigh_step = i
                break
        #MEmodel = xs.Plot.model(2)
        ME1_all = xs.Plot.addComp(1,2)
        ME2_all = xs.Plot.addComp(2,2)
        MEenergy = MEenergy_all[MElow_step:MEhigh_step]
        ME1 = ME1_all[MElow_step:MEhigh_step]
        ME2 = ME2_all[MElow_step:MEhigh_step]
        
        HEenergy_all = xs.Plot.x(3)
        for i in range(len(HEenergy_all)):
            if HEenergy_all[i] >= 27:
                HElow_step = i
                break
        #HEmodel = xs.Plot.model(3)
        HE1_all = xs.Plot.addComp(1,3)
        HE2_all = xs.Plot.addComp(2,3)
        HEenergy = HEenergy_all[HElow_step:]
        HE1 = HE1_all[HElow_step:]
        HE2 = HE2_all[HElow_step:]
        
        ax1.set_ylim(bottom=np.min(HEkeVsq)/10.)
        ax1.plot(LEenergy, np.array(LE1), 'C0', linestyle='--', linewidth=0.5, label="phabs*pow*fdcut*cyclabs*cyclabs")
        ax1.plot(MEenergy, np.array(ME1)/temp_const_m2, 'C1', linestyle='--', linewidth=0.5)
        ax1.plot(HEenergy, np.array(HE1)/temp_const_m3, 'C2', linestyle='--', linewidth=0.5)
        ax1.plot(LEenergy, np.array(LE2)/temp_phabs_nH, 'C0', linestyle='-.', linewidth=0.5, label="gauss")
        ax1.plot(MEenergy, np.array(ME2)/temp_const_m2/temp_phabs_nH, 'C1', linestyle='-.', linewidth=0.5)
        ax1.plot(HEenergy, np.array(HE2)/temp_const_m3/temp_phabs_nH, 'C2', linestyle='-.', linewidth=0.5)
        #ax1.legend(loc='upper right', fontsize=15)

        xs.Plot("del")
        LEenergy = xs.Plot.x(1)
        LEdelchi = xs.Plot.y(1)
        LExerror = xs.Plot.xErr(1)
        LEyerror = xs.Plot.yErr(1)
        
        MEenergy = xs.Plot.x(2)
        MEdelchi = xs.Plot.y(2)
        MExerror = xs.Plot.xErr(2)
        MEyerror = xs.Plot.yErr(2)
        
        HEenergy = xs.Plot.x(3)
        HEdelchi = xs.Plot.y(3)
        HExerror = xs.Plot.xErr(3)
        HEyerror = xs.Plot.yErr(3)
        
        ax2.set_xscale('log')
        #ax2.set_ylabel('(data-model)/error',size=15)
        #ax2.set_xlabel('Energy (keV)',size=20)
        ax2.set_xlim(2,105)
        ax2.set_ylim(-3,3)
        ax2.set_yticks([-2,-1,1,2])
        ax2.tick_params(axis='x', direction='in', which='both', top=True, labelbottom=False)
        ax2.tick_params(axis='y', labelsize=11, direction='out', which='both', right=True)
        
        ax2.errorbar(LEenergy, LEdelchi, yerr=LEyerror, xerr=LExerror, fmt='+', markersize=0.75, elinewidth=0.65)
        ax2.errorbar(MEenergy, MEdelchi, yerr=MEyerror, xerr=MExerror, fmt='+', markersize=0.75, elinewidth=0.65)
        ax2.errorbar(HEenergy, HEdelchi, yerr=HEyerror, xerr=HExerror, fmt='+', markersize=0.75, elinewidth=0.65)
        ax2.plot([2,105],[0,0],'g')
        ax2.text(2.1, 2.1, '(a) constant*phabs*(pow*fdcut*cyclabs*cyclabs + gaussian)', ha='left',fontsize=12)
        
        print("**************************model2**********************************") 
        xs.AllModels.clear()
        model_string = "constant*phabs*(pow*fdcut*cyclabs*cyclabs)"
        m1 = xs.Model("%s"%model_string)
        m2 = xs.AllModels(2)
        m3 = xs.AllModels(3)
        m1.setPars("%s,-1"%temp_const_m1,"%s,-1"%temp_phabs_nH,"%s,-1"%temp_pl_gamma,"%s,-1"%temp_pl_norm,"%s,-1"%temp_fdcut_cutoffE,"%s,-1"%temp_fdcut_foldE,{7:"%s,-1"%temp_cyclabs_D2},{8:"%s,-1"%temp_cyclabs_E2},{9:"%s,-1"%temp_cyclabs_W2},{12:"%s,-1"%temp_cyclabs_D1},{13:"%s,-1"%temp_cyclabs_E1},{14:"%s,-1"%temp_cyclabs_W1})
        m2.setPars("%s,-1"%temp_const_m2)#,,0.9,0.9,1.1,1.1")
        m3.setPars("%s,-1"%temp_const_m3)#,,0.9,0.9,1.1,1.1")
        xs.AllModels.show()

        xs.Plot("del")
        LEenergy = xs.Plot.x(1)
        LEdelchi = xs.Plot.y(1)
        LExerror = xs.Plot.xErr(1)
        LEyerror = xs.Plot.yErr(1)
        
        MEenergy = xs.Plot.x(2)
        MEdelchi = xs.Plot.y(2)
        MExerror = xs.Plot.xErr(2)
        MEyerror = xs.Plot.yErr(2)
        
        HEenergy = xs.Plot.x(3)
        HEdelchi = xs.Plot.y(3)
        HExerror = xs.Plot.xErr(3)
        HEyerror = xs.Plot.yErr(3)
        
        ax3.set_xscale('log')
        #ax3.set_ylabel('(data-model)/error',loc='top', size=15)
        #ax3.set_xlabel('Energy (keV)',size=20)
        ax3.set_xlim(2,105)
        ax3.set_ylim(-5.2,5.2)
        ax3.set_yticks([-4,-2,2,4])
        ax3.tick_params(axis='x', direction='in', which='both', top=True, labelbottom=False)
        ax3.tick_params(axis='y', labelsize=11, direction='out', which='both', right=True)
        
        ax3.errorbar(LEenergy, LEdelchi, yerr=LEyerror, xerr=LExerror, fmt='+', markersize=0.75, elinewidth=0.65)
        ax3.errorbar(MEenergy, MEdelchi, yerr=MEyerror, xerr=MExerror, fmt='+', markersize=0.75, elinewidth=0.65)
        ax3.errorbar(HEenergy, HEdelchi, yerr=HEyerror, xerr=HExerror, fmt='+', markersize=0.75, elinewidth=0.65)
        ax3.plot([2,105],[0,0],'g')
        ax3.text(2.1, 2.7, '(b) constant*phabs*pow*fdcut*cyclabs*cyclabs', ha='left', fontsize=12)


        print("**************************model3**********************************") 
        xs.AllModels.clear()
        model_string = "constant*phabs*(pow*fdcut*cyclabs)"
        m1 = xs.Model("%s"%model_string)
        m2 = xs.AllModels(2)
        m3 = xs.AllModels(3)
        m1.setPars("%s,-1"%temp_const_m1,"%s,-1"%temp_phabs_nH,"%s,-1"%temp_pl_gamma,"%s,-1"%temp_pl_norm,"%s,-1"%temp_fdcut_cutoffE,"%s,-1"%temp_fdcut_foldE,{7:"%s,-1"%temp_cyclabs_D2},{8:"%s,-1"%temp_cyclabs_E2},{9:"%s,-1"%temp_cyclabs_W2})
        m2.setPars("%s,-1"%temp_const_m2)#,,0.9,0.9,1.1,1.1")
        m3.setPars("%s,-1"%temp_const_m3)#,,0.9,0.9,1.1,1.1")
        xs.AllModels.show()

        xs.Plot("del")
        LEenergy = xs.Plot.x(1)
        LEdelchi = xs.Plot.y(1)
        LExerror = xs.Plot.xErr(1)
        LEyerror = xs.Plot.yErr(1)
        
        MEenergy = xs.Plot.x(2)
        MEdelchi = xs.Plot.y(2)
        MExerror = xs.Plot.xErr(2)
        MEyerror = xs.Plot.yErr(2)
        
        HEenergy = xs.Plot.x(3)
        HEdelchi = xs.Plot.y(3)
        HExerror = xs.Plot.xErr(3)
        HEyerror = xs.Plot.yErr(3)
        
        ax4.set_xscale('log')
        ax4.set_ylabel('(data-model)/error',loc='center', size=15)
        ax4.set_xlabel('Energy (keV)',size=20)
        ax4.set_xlim(2,105)
        ax4.set_yticks([-6,-2,2])
        ax4.tick_params(axis='x', direction='in', which='both', top=True, labelbottom=False)
        ax4.tick_params(axis='y', labelsize=11, direction='out', which='both', right=True)
        
        ax4.errorbar(LEenergy, LEdelchi, yerr=LEyerror, xerr=LExerror, fmt='+', markersize=0.75, elinewidth=0.65)
        ax4.errorbar(MEenergy, MEdelchi, yerr=MEyerror, xerr=MExerror, fmt='+', markersize=0.75, elinewidth=0.65)
        ax4.errorbar(HEenergy, HEdelchi, yerr=HEyerror, xerr=HExerror, fmt='+', markersize=0.75, elinewidth=0.65)
        ax4.plot([2,105],[0,0],'g')
        ax4.text(2.1, 2.7, '(c) constant*phabs*pow*fdcut*cyclabs', ha='left', fontsize=12)

        print("**************************model4**********************************") 
        xs.AllModels.clear()
        model_string = "constant*phabs*(pow*fdcut)"
        m1 = xs.Model("%s"%model_string)
        m2 = xs.AllModels(2)
        m3 = xs.AllModels(3)
        m1.setPars("%s,-1"%temp_const_m1,"%s,-1"%temp_phabs_nH,"%s,-1"%temp_pl_gamma,"%s,-1"%temp_pl_norm,"%s,-1"%temp_fdcut_cutoffE,"%s,-1"%temp_fdcut_foldE)
        m2.setPars("%s,-1"%temp_const_m2)#,,0.9,0.9,1.1,1.1")
        m3.setPars("%s,-1"%temp_const_m3)#,,0.9,0.9,1.1,1.1")
        xs.AllModels.show()

        xs.Plot("del")
        LEenergy = xs.Plot.x(1)
        LEdelchi = xs.Plot.y(1)
        LExerror = xs.Plot.xErr(1)
        LEyerror = xs.Plot.yErr(1)
        
        MEenergy = xs.Plot.x(2)
        MEdelchi = xs.Plot.y(2)
        MExerror = xs.Plot.xErr(2)
        MEyerror = xs.Plot.yErr(2)
        
        HEenergy = xs.Plot.x(3)
        HEdelchi = xs.Plot.y(3)
        HExerror = xs.Plot.xErr(3)
        HEyerror = xs.Plot.yErr(3)
        
        ax5.set_xscale('log')
        #ax5.set_ylabel('(data-model)/error',loc='center', size=15)
        ax5.set_xlabel('Energy (keV)',size=15)
        ax5.set_xlim(2,105)
        ax5.set_xticks([5,10,20,50,100])
        ax5.set_yticks([-10,-2,2])
        ax5.set_xticklabels(('5', '10', '20', '50', '100'))
        ax5.tick_params(axis='x', direction='in', which='both', top=True)
        ax5.tick_params(axis='y', labelsize=11, direction='out', which='both', right=True)
        
        ax5.errorbar(LEenergy, LEdelchi, yerr=LEyerror, xerr=LExerror, fmt='+', markersize=0.75, elinewidth=0.65)
        ax5.errorbar(MEenergy, MEdelchi, yerr=MEyerror, xerr=MExerror, fmt='+', markersize=0.75, elinewidth=0.65)
        ax5.errorbar(HEenergy, HEdelchi, yerr=HEyerror, xerr=HExerror, fmt='+', markersize=0.75, elinewidth=0.65)
        ax5.plot([2,105],[0,0],'g')
        ax5.text(2.1, 2.7, '(d) constant*phabs*pow*fdcut', ha='left', fontsize=12)

        #if not os.path.exists(path_out):
        #    os.mkdir(path_out)
        #plt.savefig('%s/%s%s.png'%(path_out,root[-3:-1],root[-1]),dpi=300)
        plt.savefig('spec_fit.eps',dpi=300)
        #plt.show()
        xs.AllModels.clear()
        xs.AllData.clear()

np.savetxt("%s/params.csv"%path_out, model_par, delimiter=",", fmt='%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s')
