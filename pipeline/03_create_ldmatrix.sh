#!/bin/bash

source PATHS
source CONFIG

RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
LD_JOBSUBDIR="${JOBSUBDIR}/ldjobs"

THIS_CONFIGDIR="${BASEDIR}/selectloci/${CONFIGDIR}"
LOCIDEFDIR="${THIS_CONFIGDIR}/locidef"
LOCUSNAMES="${THIS_CONFIGDIR}/LOCUSNAMES"
LDBASEDIR="${THIS_CONFIGDIR}/ldmap"
LDMAPWGHTDIR="${LDBASEDIR}/weighted"
LDMAPCOMBDIR="${LDBASEDIR}/combined"

for file in `ls -1 ${LOCIDEFDIR}/*.selected`; do echo `basename $file .selected`; done > ${LOCUSNAMES}


if [ ! -d ${LD_JOBSUBDIR} ];   then mkdir -p ${LD_JOBSUBDIR};   fi
cd ${LD_JOBSUBDIR}

CREATEBGEN_JOBNAME="create_bgen_${RANDSTRING}"
for STUDY in ${STUDYNAMES[@]}; do
    GENODIR=${THIS_CONFIGDIR}/${STUDY}
    JOBNAME="${CREATEBGEN_JOBNAME}_${STUDY}"
    sed "s|_JOBNAME|${JOBNAME}|g;
         s|_GENDIR_|${GENODIR}|g;
         s|_LOC_NM_|${LOCUSNAMES}|g;
         s|_QC_TOL_|${QCTOOL}|g;
        " ${MASTER_BSUBDIR}/createbgen.bsub > ${JOBNAME}.bsub
    bsub < ${JOBNAME}.bsub
done

LDSTORE_JOBNAME="ldmatrix_${RANDSTRING}"
for STUDY in ${STUDYNAMES[@]}; do
    JOBNAME="${LDSTORE_JOBNAME}_${STUDY}"
    OUTDIR="${LDBASEDIR}/${STUDY}"
    sed "s|_JOBNAME|${JOBNAME}|g;
         s|_GSTUDY_|${STUDY}|g;
         s|_USELOCI|${LOCUSNAMES}|g;
         s|_LOCIDIR|${THIS_CONFIGDIR}/${STUDY}|g;
         s|_LDSTORE|${LDSTORE}|g;
         s|_OUTDIR_|${OUTDIR}|g;
         " ${MASTER_BSUBDIR}/ldstore.bsub > ${JOBNAME}.bsub
    bsub -w "done(${CREATEBGEN_JOBNAME}*)" < ${JOBNAME}.bsub
done

while read LOCUSPREFIX; do
    WGT_LD_JOBNAME="weighted_LD_${LOCUSPREFIX}"
    sed -e "s|_JOBNAME|${WGT_LD_JOBNAME}|g;
            s|_WGHT_LD|${LDMAP_WEIGHTED}|g;
            s|_LOCUSP_|${LOCUSPREFIX}|g;
            s|_STUDYN_|\"${STUDYNAMES[*]}\"|g;
            s|_SAMPLES|\"${STUDYSAMPLES[*]}\"|g;
            s|_METDIR_|${LOCIDEFDIR}|g;
            s|_LDBASE_|${LDBASEDIR}|g;
            s|_LD_DIR_|${LDMAPWGHTDIR}|g;
            " ${MASTER_BSUBDIR}/ldmap_weighted.bsub > ${WGT_LD_JOBNAME}.bsub
    bsub -w "done(${LDSTORE_JOBNAME}*)" < ${WGT_LD_JOBNAME}.bsub
done < ${LOCUSNAMES}

cd ${CURDIR}
