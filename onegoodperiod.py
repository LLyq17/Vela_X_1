from scipy.optimize import leastsq
import numpy as np
from astropy.io import fits
import os
energy=["HE.fes","ME.fes","LE.fes"]
def gau(x,mu,sigma):
    result=1./np.sqrt(2*np.pi)/sigma*np.exp(-(x-mu)**2./(2.*sigma**2.))
    return result
def f_err(p,y,x):
    return y-gau(x,*p)
#result=[]
#error=[]
#time=[]
os.system('rm ./out/period_* && rm period_* && rm freq_*')
MJDREFF = 0.00076601852000000
MJDREFI = 55927

ObsID=["P020101214106-20191102-01-01", "P020101214107-20191103-02-01", "P020101214108-20191103-02-01", "P020101214209-20191121-02-01", "P020101214210-20191121-02-01", "P020101214301-20191215-01-01", "P020101214302-20191216-02-01", "P020101214401-20191216-01-01", "P020101214402-20191216-01-01", "P020101214403-20191217-02-01", "P020101214405-20191217-02-01", "P020101214501-20200103-01-01", "P020101214503-20200103-01-01", "P020101214504-20200103-01-01", "P020101214510-20200104-02-01", "P020101214704-20200122-01-01", "P020101214705-20200122-01-01", "P020101214707-20200122-01-01", "P020101214708-20200123-02-01", "P020101214709-20200123-02-01", "P020101214711-20200123-02-01", "P020101214803-20200207-01-01", "P020101214901-20200322-01-01", "P020101214903-20200322-01-01", "P020101214906-20200323-02-01", "P020101214910-20200323-02-01", "P020101214911-20200323-02-01", "P020101214913-20200323-02-01", "P020101214914-20200324-03-01", "P020101214915-20200324-03-01", "P020101215002-20200421-02-01", "P020101215103-20200617-01-01", "P020101215203-20200623-01-01"]

for i in range(len(ObsID)):
    result=[]
    error=[]
    time=[]
    result_f=[]
    error_f=[]
    filelist="/home/lq/hxmt/Vela_X-1/P0201012_v204/Period/%s/"%(ObsID[i])
    for e in range(2):
        try:
          newfile=filelist+"/"+energy[e]
          data=fits.open(newfile)
        except:
          print("error occur!"+"job"+str(i))
        if e == 0:
          chisq=data[1].data.field("CHISQRD1")
        else:
          chisq=chisq+data[1].data.field("CHISQRD1")
    p=data[1].data.field("PERIOD")
    c, ret_val=leastsq(f_err,[283.,0.1],args=(chisq,p))
    result.append(c[0])
    result_f.append(1/c[0])
    time.append(data[1].header["TSTARTI"]+40000+data[1].header["TSTARTF"])  #mjd
    #time.append((data[1].header["TSTARTI"]+40000+data[1].header["TSTARTF"]- MJDREFI - MJDREFF)*86400)
    error.append(c[1])
    error_f.append(1/c[0]-1/(c[1]+c[0]))
    data.close()

    xy=np.array([time,result]).T
    xy_f=np.array([time,result_f]).T
    np.savetxt("./out/period_%s.txt"%(ObsID[i].split("-")[0]),xy)
    os.system('cat ./out/period_%s.txt >>period_value.txt'%(ObsID[i].split("-")[0]))
    np.savetxt("./out/period_err%s.txt"%(ObsID[i].split("-")[0]),error)
    os.system('cat ./out/period_err%s.txt >>period_err.txt'%(ObsID[i].split("-")[0]))

    np.savetxt("./out/freq_%s.txt"%(ObsID[i].split("-")[0]),xy_f)
    os.system('cat ./out/freq_%s.txt >>freq_value.txt'%(ObsID[i].split("-")[0]))
    np.savetxt("./out/freq_err%s.txt"%(ObsID[i].split("-")[0]),error_f)
    os.system('cat ./out/freq_err%s.txt >>freq_error.txt'%(ObsID[i].split("-")[0]))

