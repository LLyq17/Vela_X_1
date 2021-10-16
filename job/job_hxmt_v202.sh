#!/bin/bash
#SBATCH -J hxmt
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --mem-per-cpu=5GB
#SBATCH -p batch
#SBATCH --output="LOG/job_JOBID_%j.log"
#SBATCH --error="ERR/job_JOBID_%j.err"

# Bourne shell (sh/bash):
export HEADASNOQUERY=
export HEADASPROMPT=/dev/null

export PATH=/data/Python2/bin/:$PATH && export HEADAS=/data/software/hxmtsoftv2.02_install/x86_64-pc-linux-gnu-libc2.17 && . $HEADAS/headas-init.sh && export CALDBCONFIG=/data/software/hxmtsoftv2.02_install/CALDB2.02/caldb.config && export CALDB=/data/software/hxmtsoftv2.02_install/CALDB2.02

# 低能LE数据预处理
lepical evtfile=FULL/LE/HXMT_PART_LE-Evt_FFFFFF_V1_L1P.FITS gainfile=CALDB maxtimedel=60 \
tempfile=FULL/LE/HXMT_PART_LE-TH_FFFFFF_V1_L1P.FITS  outfile=OUT/lepi.fits

lerecon evtfile=OUT/lepi.fits outfile=OUT/lerecon.fits instatusfile=FULL/LE/HXMT_PART_LE-InsStat_FFFFFF_V1_L1P.FITS

legtigen evtfile=NONE instatusfile=FULL/LE/HXMT_PART_LE-InsStat_FFFFFF_V1_L1P.FITS \
tempfile=FULL/LE/HXMT_PART_LE-TH_FFFFFF_V1_L1P.FITS ehkfile=FULL/AUX/HXMT_PART_EHK_FFFFFF_V1_L1P.FITS \
outfile=OUT/legti.fits defaultexpr=NONE expr='ELV>10&&COR>8&&DYE_ELV>30&&SAA_FLAG==0&&ANG_DIST<0.04&&T_SAA>300&&TN_SAA>300'

legti OUT/lerecon.fits OUT/legti.fits OUT/lenewgti.fits

lescreen evtfile=OUT/lerecon.fits gtifile=OUT/lenewgti.fits outfile=OUT/lescreen_all.fits userdetid="0-94"

hxbary evtfile=OUT/lescreen_all.fits orbitfile=FULL/ACS/HXMT_PART_Orbit_FFFFFF_V1_L1P.FITS ra=RA dec=DEC eph=2

lelcgen minPI=106 maxPI=1170 evtfile=OUT/lescreen_all.fits outfile=OUT/lelc binsize=0.0078125 \
userdetid="0,2-4,6-10,12,14,20,22-26,28-30,32,34-36,38-42,44,46,52,54-58,60-62,64,66-68,70-74,76,78,84,86-90,92-94" eventtype=1

ls OUT/lelc_g0_0-94.lc > OUT/lelc.txt

lebkgmap lc OUT/lescreen_all.fits OUT/lenewgti.fits OUT/lelc.txt 106 1170 OUT/lelcbkg

lespecgen evtfile=OUT/lescreen_all.fits outfile=OUT/lespec \
userdetid="0,2-4,6-10,12,14,20,22-26,28,30,32,34-36,38-42,44,46,52,54-58,60-62,64,66-68,70-74,76,78,84,86-90,92-94" eventtype=1

ls OUT/lespec_g0_0-94.pha > OUT/lespec.txt

lebkgmap spec OUT/lescreen_all.fits OUT/lenewgti.fits OUT/lespec.txt 0 1535 OUT/lespecbkg

lerspgen phafile=OUT/lespec_g0_0-94.pha attfile=FULL/ACS/HXMT_PART_Att_FFFFFF_V1_L1P.FITS \
tempfile=FULL/LE/HXMT_PART_LE-TH_FFFFFF_V1_L1P.FITS ra=RA dec=DEC outfile=OUT/lersp.fits

fparkey lespecbkg.pha OUT/lespec_g0_0-94.pha BACKFILE 
fparkey lersp.fits OUT/lespec_g0_0-94.pha RESPFILE

