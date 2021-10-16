#!/bin/bash
#SBATCH -J hxmt
#SBATCH -N 1 
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=5GB
#SBATCH -p batch
#SBATCH --output="LOG/job_JOBID_%j.log"
#SBATCH --error="ERR/job_JOBID_%j.err"

#export OMP_NUM_THREADS=6

# Bourne shell (sh/bash):
export HEADASNOQUERY=
export HEADASPROMPT=/dev/null

export PATH=/data/Python2/bin/:$PATH && export HEADAS=/data/software/hxmtsoftv2.02_install/x86_64-pc-linux-gnu-libc2.17 && . $HEADAS/headas-init.sh && export CALDBCONFIG=/data/software/hxmtsoftv2.02_install/CALDB2.02/caldb.config && export CALDB=/data/software/hxmtsoftv2.02_install/CALDB2.02

# 高能HE数据预处理**********************************************************************
hepical evtfile=`ls FULL/HE/HXMT_PART_HE-Evt_FFFFFF_V[1-9]_L1P.FITS` outfile=OUT/hepi.fits

hegtigen hvfile=`ls FULL/HE/HXMT_PART_HE-HV_FFFFFF_V[1-9]_L1P.FITS` tempfile=`ls FULL/HE/HXMT_PART_HE-TH_FFFFFF_V[1-9]_L1P.FITS` ehkfile=`ls FULL/AUX/HXMT_PART_EHK_FFFFFF_V[1-9]_L1P.FITS` outfile=OUT/hegti.fits defaultexpr=NONE expr='ELV>10&&COR>8&&SAA_FLAG==0&&ANG_DIST<0.04&&T_SAA>300&&TN_SAA>300' pmfile=`ls FULL/HE/HXMT_PART_HE-PM_FFFFFF_V[1-9]_L1P.FITS`

hescreen evtfile=OUT/hepi.fits gtifile=OUT/hegti.fits outfile=OUT/hescreen_all.fits userdetid="0-17" eventtype=1 anticoincidence=yes starttime=0 stoptime=0 minPI=0 maxPI=255

hxbary evtfile=OUT/hescreen_all.fits orbitfile=`ls FULL/ACS/HXMT_PART_Orbit_FFFFFF_V[1-9]_L1P.FITS` ra=RA dec=DEC eph=2

helcgen deadcorr=yes minPI=26 maxPI=150 evtfile=OUT/hescreen_all.fits deadfile=`ls FULL/HE/HXMT_PART_HE-DTime_FFFFFF_V[1-9]_L1P.FITS` outfile=OUT/helc binsize=0.0078125 userdetid="0-15,17"

ls OUT/helc_g0_0-17.lc > OUT/helc.txt

hebkgmap lc OUT/hescreen_all.fits FULL/AUX/HXMT_PART_EHK_FFFFFF_V[1-9]_L1P.FITS OUT/hegti.fits FULL/HE/HXMT_PART_HE-DTime_FFFFFF_V[1-9]_L1P.FITS OUT/helc.txt 26 150 OUT/helcbkg

hespecgen evtfile=OUT/hescreen_all.fits deadfile=`ls FULL/HE/HXMT_PART_HE-DTime_FFFFFF_V[1-9]_L1P.FITS` outfile=OUT/hespec userdetid="0;1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17" eventtype=1 starttime=0 stoptime=0 minPI=0 maxPI=255

ls OUT/hespec_g* | sort -V > OUT/hespec.txt

hebkgmap spec OUT/hescreen_all.fits FULL/AUX/HXMT_PART_EHK_FFFFFF_V[1-9]_L1P.FITS OUT/hegti.fits FULL/HE/HXMT_PART_HE-DTime_FFFFFF_V[1-9]_L1P.FITS OUT/hespec.txt 0 255 OUT/hespecbkg

herspgen phafile=OUT/hespec_g17_17.pha attfile=`ls FULL/ACS/HXMT_PART_Att_FFFFFF_V[1-9]_L1P.FITS` ra=RA dec=DEC outfile=OUT/hersp_g17_17.fits

fparkey hespecbkg_17.pha OUT/hespec_g17_17.pha BACKFILE
fparkey hersp_g17_17.fits OUT/hespec_g17_17.pha RESPFILE

for ((i=0;i<=15;i++));
do

herspgen phafile=OUT/hespec_g$i\_$i\.pha attfile=`ls FULL/ACS/HXMT_PART_Att_FFFFFF_V[1-9]_L1P.FITS` ra=RA dec=DEC outfile=OUT/hersp_g$i\_$i\.fits

fparkey hespecbkg_$i.pha OUT/hespec_g$i\_$i\.pha BACKFILE
fparkey hersp_g$i\_$i\.fits OUT/hespec_g$i\_$i\.pha RESPFILE

done

ls OUT/hespec_g* | sort -V > OUT/src.dat
ls OUT/hespecbkg_* | sort -V > OUT/bkg.dat
ls OUT/hersp_g* | sort -V > OUT/rsp.dat
sed -i '/hespec_g16_16.pha/d' OUT/src.dat
sed -i '/hespecbkg_16.pha/d' OUT/bkg.dat
hhe_spec2pi OUT/src.dat OUT/bkg.dat OUT/rsp.dat OUT/out_spec_he.pha OUT/outbgd_he.pha OUT/out_response_he.pha
fparkey outbgd_he.pha OUT/out_spec_he.pha BACKFILE
fparkey out_response_he.pha OUT/out_spec_he.pha RESPFILE


# 中能ME数据预处理***************************************************************************
mepical evtfile=`ls FULL/ME/HXMT_PART_ME-Evt_FFFFFF_V[1-9]_L1P.FITS` tempfile=`ls FULL/ME/HXMT_PART_ME-TH_FFFFFF_V[1-9]_L1P.FITS` outfile=OUT/mepi.fits

megrade evtfile=OUT/mepi.fits outfile=OUT/megrade.fits deadfile=OUT/medead.fits binsize=0.0078125 clobber=yes history=yes

megtigen tempfile=`ls FULL/ME/HXMT_PART_ME-TH_FFFFFF_V[1-9]_L1P.FITS` ehkfile=`ls FULL/AUX/HXMT_PART_EHK_FFFFFF_V[1-9]_L1P.FITS` outfile=OUT/megti.fits defaultexpr=NONE expr='ELV>10&&COR>8&&SAA_FLAG==0&&ANG_DIST<0.04&&T_SAA>300&&TN_SAA>300'

megti OUT/megrade.fits OUT/megti.fits OUT/menewgti.fits $HEADAS/refdata/medetectorstatus.fits OUT/out_detstatus.fits ##增加的参数是输入坏探测器文件，并根据观测时间产生新的对应的坏探测器文件，供后续使用

mescreen evtfile=OUT/megrade.fits gtifile=OUT/menewgti.fits outfile=OUT/mescreen_all.fits userdetid="0-53" baddetfile=OUT/out_detstatus.fits

hxbary evtfile=OUT/mescreen_all.fits orbitfile=`ls FULL/ACS/HXMT_PART_Orbit_FFFFFF_V[1-9]_L1P.FITS` ra=RA dec=DEC eph=2

melcgen deadcorr=yes minPI=119 maxPI=546 evtfile=OUT/mescreen_all.fits deadfile=OUT/medead.fits outfile=OUT/melc binsize=0.0078125 userdetid="0-7,11-17,18-25,29-35,36-43,47-53"

ls OUT/melc_g0_0-53.lc > OUT/melc.txt

mebkgmap lc OUT/mescreen_all.fits FULL/AUX/HXMT_PART_EHK_FFFFFF_V[1-9]_L1P.FITS OUT/menewgti.fits OUT/medead.fits FULL/ME/HXMT_PART_ME-TH_FFFFFF_V[1-9]_L1P.FITS OUT/melc.txt 119 546 OUT/melcbkg OUT/out_detstatus.fits

mespecgen evtfile=OUT/mescreen_all.fits deadfile=OUT/medead.fits outfile=OUT/mespec userdetid="0-7,11-25,29-43,47-53"

ls OUT/mespec_g0_0-53.pha > OUT/mespec.txt

mebkgmap spec OUT/mescreen_all.fits FULL/AUX/HXMT_PART_EHK_FFFFFF_V[1-9]_L1P.FITS OUT/menewgti.fits OUT/medead.fits FULL/ME/HXMT_PART_ME-TH_FFFFFF_V[1-9]_L1P.FITS OUT/mespec.txt 0 1023 OUT/mespecbkg OUT/out_detstatus.fits

merspgen phafile=OUT/mespec_g0_0-53.pha attfile=`ls FULL/ACS/HXMT_PART_Att_FFFFFF_V[1-9]_L1P.FITS` ra=RA dec=DEC outfile=OUT/mersp_g0_53.fits

fparkey mespecbkg.pha OUT/mespec_g0_0-53.pha BACKFILE
fparkey mersp_g0_53.fits OUT/mespec_g0_0-53.pha RESPFILE


# 低能LE数据预处理***************************************************************************
lepical evtfile=`ls FULL/LE/HXMT_PART_LE-Evt_FFFFFF_V[1-9]_L1P.FITS` gainfile=CALDB maxtimedel=60 \
tempfile=`ls FULL/LE/HXMT_PART_LE-TH_FFFFFF_V[1-9]_L1P.FITS`  outfile=OUT/lepi.fits

lerecon evtfile=OUT/lepi.fits outfile=OUT/lerecon.fits instatusfile=`ls FULL/LE/HXMT_PART_LE-InsStat_FFFFFF_V[1-9]_L1P.FITS`

legtigen evtfile=NONE instatusfile=`ls FULL/LE/HXMT_PART_LE-InsStat_FFFFFF_V[1-9]_L1P.FITS` \
tempfile=`ls FULL/LE/HXMT_PART_LE-TH_FFFFFF_V[1-9]_L1P.FITS` ehkfile=`ls FULL/AUX/HXMT_PART_EHK_FFFFFF_V[1-9]_L1P.FITS` \
outfile=OUT/legti.fits defaultexpr=NONE expr='ELV>10&&COR>8&&DYE_ELV>30&&SAA_FLAG==0&&ANG_DIST<0.04&&T_SAA>300&&TN_SAA>300'

legti OUT/lerecon.fits OUT/legti.fits OUT/lenewgti.fits

lescreen evtfile=OUT/lerecon.fits gtifile=OUT/lenewgti.fits outfile=OUT/lescreen_all.fits userdetid="0-95"

hxbary evtfile=OUT/lescreen_all.fits orbitfile=`ls FULL/ACS/HXMT_PART_Orbit_FFFFFF_V[1-9]_L1P.FITS` ra=RA dec=DEC eph=2

lelcgen minPI=106 maxPI=1170 evtfile=OUT/lescreen_all.fits outfile=OUT/lelc binsize=0.0078125 \
userdetid="0,2-4,6-10,12,14,20,22-26,28-30,32,34-36,38-42,44,46,52,54-58,60-62,64,66-68,70-74,76,78,84,86,88-90,92-94" eventtype=1

ls OUT/lelc_g0_0-94.lc > OUT/lelc.txt

lebkgmap lc OUT/lescreen_all.fits OUT/lenewgti.fits OUT/lelc.txt 106 1170 OUT/lelcbkg

lespecgen evtfile=OUT/lescreen_all.fits outfile=OUT/lespec \
userdetid="0,2-4,6-10,12,14,20,22-26,28,30,32,34-36,38-42,44,46,52,54-58,60-62,64,66-68,70-74,76,78,84,86,88-90,92-94" eventtype=1

ls OUT/lespec_g0_0-94.pha > OUT/lespec.txt

lebkgmap spec OUT/lescreen_all.fits OUT/lenewgti.fits OUT/lespec.txt 0 1535 OUT/lespecbkg

lerspgen phafile=OUT/lespec_g0_0-94.pha attfile=`ls FULL/ACS/HXMT_PART_Att_FFFFFF_V[1-9]_L1P.FITS` \
tempfile=`ls FULL/LE/HXMT_PART_LE-TH_FFFFFF_V[1-9]_L1P.FITS` ra=RA dec=DEC outfile=OUT/lersp.fits

fparkey lespecbkg.pha OUT/lespec_g0_0-94.pha BACKFILE 
fparkey lersp.fits OUT/lespec_g0_0-94.pha RESPFILE

