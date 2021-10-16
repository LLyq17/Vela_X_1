from xspec import *
import matplotlib.pyplot as plt
import math
import numpy as np
import os
import astropy.io.fits as fits
from astropy.io.fits import getheader

AllModels.initpackage("mymodel","lmodel.dat","/home/lq/heasoft-6.28/Xspec/src/Qifunc")
AllModels.lmod("mymodel","/home/lq/heasoft-6.28/Xspec/src/Qifunc")

AllData.clear()
AllData("1:1 lespec_g0_0-94.pha 2:2 mespec_g0_0-53.pha 3:3 hespec_g0_0-17.pha") 
s1 = AllData(1)
s2 = AllData(2)
s3 = AllData(3)

s1.ignore("**-2. 10.-**")
s2.ignore("**-9. 29.-**")
s3.ignore("**-27. 105.-**")

AllModels.clear()
m1 = Model("phabs*pow*fdcut*cyclabs*cyclabs*constant")
m2 = AllModels(2)
m3 = AllModels(3)

# 1:nH(phabs)  2:PhoIndex  3:norm  4:cutoffE  5:foldE
# 6:Depth0   7:E0   8:Width0    9:Depth2  10:Width2  (cyclabs_2)
# 11:Depth0  12:E0  13:Width0  14:Depth2  15:Width2   (cyclabs_1)
par7 = m1(7)
par7.values = [45., .01, 35, 35, 55, 55]
par8 = m1(8)
par8.values = 9.3
#par6 = m1(6)
#par6.values = [1.3, .01, 0, 0, 10,10]
m2(16).link = ""
m3(16).link = ""
m1(16).frozen = True
m2(16).values = [1., .001, .85, .85, 1.05, 1.05]
m3(16).values = [1., .001, .85, .85, 1.05, 1.05]

par12= m1(12)
par12.values = [25, .01, 18, 18, 30, 30]
par13= m1(13)
par13.values = 3.9
#par11 = m1(11)
#par11.values = [0.3, .01, 0, 0, 10,10]
AllModels.show()

Fit.query = "yes"
Fit.perform()
Fit.query = "yes"
Fit.nIterations = 1000000
Fit.query = "yes"
Fit.show()

par8.frozen = False
par13.frozen = False
AllModels.show()

Fit.query = "yes"
Fit.perform()
Fit.query = "yes"
Fit.nIterations = 1000000
Fit.query = "yes"
Fit.show()

Plot.addCommand("label top")
Plot.addCommand("time off")
Plot.device = "01.eps/cps"
Plot.xAxis = "KeV"
Plot.setRebin(12, 12, 1)
Plot.setRebin(15, 15, 2)
Plot.setRebin(1, 1, 3)
Plot.background = False
Plot("uf del")

Plot.device = "/xw"
Plot.xAxis = "KeV"
Plot.setRebin(12, 12, 1)
Plot.setRebin(15, 15, 2)
Plot.setRebin(1, 1, 3)
Plot.background = False
Plot("ldata euf del")

#parallel
Xset.parallel.leven = 4
Xset.parallel.error = 4
Xset.parallel.walkers = 50
Xset.parallel.show()

os.system('rm chain2.fits')
AllChains.clear()
AllChains.defBurn = 200
AllChains.defLength = 100000
#AllChains.defWalkers = 50
AllChains.defProposal = "gaussian fit"
c2 = Chain("chain2.fits")
AllChains.show()

#hdul= fits.open('chain2.fits')
#hdul[1].data
#hdr = hdul[1].header
#Plot("Chain 0")
#data=hdul[1].data[100:]
#if os.path.exists('chain1.fits'): os.system('rm chain1.fits')
#fits.writeto('chain1.fits', data, hdr)

#AllChains.clear()
#AllChains += "chain1.fits"
#AllChains.show()

Fit.error("2.706 1")
Fit.error("2.706 2")
Fit.error("2.706 4")
Fit.error("2.706 5")
Fit.error("2.706 7")
Fit.error("2.706 6")
Fit.error("2.706 8")
Fit.error("2.706 12")
Fit.error("2.706 11")
Fit.error("2.706 13")

par1 = AllModels(1)(1)
par2 = AllModels(1)(2)
par4 = AllModels(1)(4)
par5 = AllModels(1)(5)
par6 = AllModels(1)(6)
par7 = AllModels(1)(7)
par8 = AllModels(1)(8)
par11 = AllModels(1)(11)
par12 = AllModels(1)(12)
par13 = AllModels(1)(13)

print(par1.error,par2.error,par4.error,par5.error,par6.error,par7.error,par8.error,par11.error,par12.error,par13.error)

cons2=m2.constant.factor.values[0]
cons3=m3.constant.factor.values[0]
nH=m1(1).values[0]
PhoIndex=m1(2).values[0]
norm=m1(3).values[0]
cutoffE=m1(4).values[0]
foldE=m1(5).values[0]
Depth0=m1(6).values[0]
E0 =m1(7).values[0]
Width0=m1(8).values[0]
Depth0_2=m1(11).values[0]
E0_2 =m1(12).values[0]
Width0_2=m1(13).values[0]
print(norm)

f = open("para2_cyclo.txt",'w')
f.write('E2'+' '+'E2err-'+' '+'E2err+'+' '+'depth2'+' '+'d2err-'+' '+'d2err+'+' '+
        'width2'+' '+'w2err-'+' '+'w2err+'+' '+
        'E1'+' '+'E1-'+' '+'E1+'+' '+'depth1'+' '+'d1-'+' '+'d1+'+' '+
        'width1'+' '+'w1-'+' '+'w1+'+'\n')
it=[7,6,8,12,11,13,]
for i in it:
    f.write(str(round(m1(i).values[0],2))+' '+str(round((m1(i).error[0]-m1(i).values[0]),2))+' '+
            str(round((m1(i).error[1]-m1(i).values[0]),2))+' ')
f.write('\n')

f2 = open("para2_conti.txt",'w')
f2.write('nH'+' '+'nHerr-'+' '+'nHerr+'+' '+'gindex'+' '+'gerr-'+' '+'gerr+'+' '+
        'Ecut'+' '+'Ecuterr-'+' '+'Ecuterr+'+' '+'Efold'+' '+'Eferr-'+' '+'Eferr+'+'\n')
it=[1,2,4,5,]
for i in it:
    f2.write(str(round(m1(i).values[0],2))+' '+str(round((m1(i).error[0]-m1(i).values[0]),2))+' '+
            str(round((m1(i).error[1]-m1(i).values[0]),2))+' ')
f2.write('\n')

f = open("para2_cyclo.txt",'w')
f.write('E2'+' '+'E2err-'+' '+'E2err+'+' '+'depth2'+' '+'d2err-'+' '+'d2err+'+' '+
        'width2'+' '+'w2err-'+' '+'w2err+'+' '+
        'E1'+' '+'E1-'+' '+'E1+'+' '+'depth1'+' '+'d1-'+' '+'d1+'+' '+
        'width1'+' '+'w1-'+' '+'w1+'+'\n')
it=[7,6,8,12,11,13,]
for i in it:
    f.write(str(round(m1(i).values[0],2))+' '+str(round((m1(i).error[0]-m1(i).values[0]),2))+' '+
            str(round((m1(i).error[1]-m1(i).values[0]),2))+' ')
f.write('\n')

f2 = open("para2_conti.txt",'w')
f2.write('nH'+' '+'nHerr-'+' '+'nHerr+'+' '+'gindex'+' '+'gerr-'+' '+'gerr+'+' '+
        'Ecut'+' '+'Ecuterr-'+' '+'Ecuterr+'+' '+'Efold'+' '+'Eferr-'+' '+'Eferr+'+'\n')
it=[1,2,4,5,]
for i in it:
    f2.write(str(round(m1(i).values[0],2))+' '+str(round((m1(i).error[0]-m1(i).values[0]),2))+' '+
            str(round((m1(i).error[1]-m1(i).values[0]),2))+' ')
f2.write('\n')

AllModels.clear()
m1 = Model("constant*phabs*cyclabs*cyclabs*cflux(pow*fdcut)")

m2 = AllModels(2)
m3 = AllModels(3)

m1.phabs.nH.values=nH
m1.phabs.nH.frozen = True
m1.powerlaw.PhoIndex.values=PhoIndex
m1.powerlaw.PhoIndex.frozen = True
m1.powerlaw.norm.values=norm
m1.powerlaw.norm.frozen = True
m1.fdcut.cutoffE.values=cutoffE
m1.fdcut.cutoffE.frozen = True
m1.fdcut.foldE.values=foldE
m1.fdcut.foldE.frozen = True
m1.cyclabs.Depth0.values=Depth0
m1.cyclabs.Depth0.frozen = True
m1.cyclabs.E0.values=E0
m1.cyclabs.E0.frozen = True
m1.cyclabs.Width0.values=Width0
m1.cyclabs.Width0.frozen = True
m1(8).values=Depth0_2
m1(8).frozen = True
m1(9).values=E0_2
m1(9).frozen = True
m1(10).values=Width0_2
m1(10).frozen = True

m1.cflux.Emin.values=2.
m1.cflux.Emax.values=10.
m1.cflux.lg10Flux.values= -8.4
m2.cflux.Emin.link = ""
m3.cflux.Emin.link = ""
m2.cflux.Emin.values=10.
m2.cflux.Emax.values=29.
m2.cflux.lg10Flux.values= -8.1
m3.cflux.Emin.values=29.
m3.cflux.Emax.values=105.
m3.cflux.lg10Flux.values= -8.5

m2.constant.factor.link = ""
m3.constant.factor.link = ""
m2.constant.factor.values=cons2
m3.constant.factor.values=cons3
m1.constant.factor.frozen = True
m2.constant.factor.frozen = True
m3.constant.factor.frozen = True

AllModels.show()

Fit.perform()
Fit.nIterations = 1000000
Fit.show()

Plot.device = "/xw"
Plot.xAxis = "KeV"
Plot.setRebin(12, 12, 1)
Plot.setRebin(15, 15, 2)
Plot.setRebin(1, 1, 3)
Plot.background = False
Plot("ldata euf del")

os.system('rm chain2.fits')
AllChains.clear()
AllChains.defBurn = 200
AllChains.defLength = 100000
#AllChains.defWalkers = 50
AllChains.defProposal = "gaussian fit"
c2 = Chain("chain2.fits")
AllChains.show()

#hdul= fits.open('chain2.fits')
#hdul[1].data
#hdr = hdul[1].header
#Plot("Chain 0")
#data=hdul[1].data[100:]
#if os.path.exists('chain1.fits'): os.system('rm chain1.fits')
#fits.writeto('chain1.fits', data, hdr)

#AllChains.clear()
#AllChains += "chain1.fits"
#AllChains.show()

Fit.error("2.706 15")
Fit.error("2.706 34")
Fit.error("2.706 53")

par15 = AllModels(1)(15)
par15.error

par15 = AllModels(1)(15)
par34 = AllModels(2)(15)
par53 = AllModels(3)(15)

print(par15.values[0],par34.values[0],par53.values[0])
print(par15.error,par34.error,par53.error)

s_10=pow(10,par15.values[0])+ pow(10,par34.values[0])+ pow(10,par53.values[0])
print(s_10)
s=np.log10(s_10)
serr_1= np.sqrt((par15.error[0]-par15.values[0])**2+(par34.error[0]-par34.values[0])**2+(par53.error[0]-par53.values[0])**2)
serr_2= np.sqrt((par15.error[1]-par15.values[0])**2+(par34.error[1]-par34.values[0])**2+(par53.error[1]-par53.values[0])**2)

f3 = open("para2_cflux.txt",'w')
f3.write('flux'+' '+'fluxerr-'+' '+'fluxerr+'+'\n')

f3.write(str(round(s,12))+' '+str(round(serr_1,12))+' '+str(round(serr_2,12))+'\n')

f3 = open("para2_cflux.txt",'w')
f3.write('flux'+' '+'fluxerr-'+' '+'fluxerr+'+'\n')

f3.write(str(round(s,12))+' '+str(round(serr_1,12))+' '+str(round(serr_2,12))+'\n')



