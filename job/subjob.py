import os, sys
# writed by liuqi at 2021.01.27

NUM=sys.argv[1]  
samplePath="/data/sci_data/hxmt/%s"%NUM
os.system("mkdir -vp %s"%NUM.split("/")[-1])
outputPath=os.path.abspath(NUM.split("/")[-1])
print("Your input path is: "+samplePath)
print("Your output path is: "+outputPath)
os.chdir(outputPath)
print(os.getcwd())

for i in os.listdir(samplePath):
    if os.path.exists(samplePath+"/"+"LE"):
        os.system("mkdir HE ME LE")
        break
    else:
        if i!="ACS" and  i!="AUX" and i!="ExpoList.XML" and i!="FileList.FITS":
            os.system("mkdir -vp %s/%s"%(outputPath,i))
            os.chdir(outputPath+"/"+i)

            for j in os.listdir(samplePath+"/"+i):
                if os.path.exists(samplePath+"/"+i+"/"+"LE"):
                    os.system("mkdir HE ME LE")
                    break
                else:
                    if j!="ACS" and  j!="AUX" and j!="ExpoList.XML" and j!="FileList.FITS":
                        os.system("mkdir -vp %s/%s/%s"%(outputPath,i,j))
                        os.chdir(outputPath+"/"+i+"/"+j)
                        #os.system("mkdir HE ME LE")
                        #os.chdir(outputPath+"/"+i+"/"+j+"/"+"LE")
                        os.chdir(outputPath)
                        os.system("mkdir -vp config/%s log/%s err/%s"%(i,i,i))
                        os.system('cp /data/users/liuqi/data/hxmt/job_hxmt_v202.sh config/%s/job_hxmt_v202_%s.sh'%(i,j))
                        os.system('sed -i "s#FULL#%s#g" config/%s/job_hxmt_v202_%s.sh'%(samplePath+"/"+i+"/"+j,i,j))
                        os.system('sed -i "s#PART#%s#g" config/%s/job_hxmt_v202_%s.sh'%(j.split("-")[0],i,j)) 
                        os.system('sed -i "s#RA#%s#g" config/%s/job_hxmt_v202_%s.sh'%(sys.argv[2],i,j))
                        os.system('sed -i "s#DEC#%s#g" config/%s/job_hxmt_v202_%s.sh'%(sys.argv[3],i,j)) 
                        os.system('sed -i "s#OUT#%s#g" config/%s/job_hxmt_v202_%s.sh'%(outputPath+"/"+i+"/"+j,i,j))                        
                        os.system('sed -i "s#LOG#log/%s#g" config/%s/job_hxmt_v202_%s.sh'%(i,i,j))
                        os.system('sed -i "s#ERR#err/%s#g" config/%s/job_hxmt_v202_%s.sh'%(i,i,j))
                        os.system('sed -i "s#JOBID#%s#g" config/%s/job_hxmt_v202_%s.sh'%(j,i,j))
                        os.system('sbatch config/%s/job_hxmt_v202_%s.sh'%(i,j))
