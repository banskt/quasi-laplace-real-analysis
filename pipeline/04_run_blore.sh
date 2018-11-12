#!/bin/bash

source PATHS
source CONFIG

RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
BLORE_JOBSUBDIR="${JOBSUBDIR}/${CONFIGDIR}/blore"

THIS_CONFIGDIR="${BASEDIR}/${CONFIGDIR}"
THIS_LOCIDIR="${THIS_CONFIGDIR}/selectloci"
SAMPLEDIR="${BASEDIR}/metaanalysis/samples"
USELOCI="${THIS_CONFIGDIR}/LOCUSNAMES"
BLORE_RESDIR="${THIS_CONFIGDIR}/blore/summary_stat"
BLORE_METARS="${THIS_CONFIGDIR}/blore/meta"

if [ ! -d ${BLORE_JOBSUBDIR} ]; then mkdir -p ${BLORE_JOBSUBDIR}; fi
cd ${BLORE_JOBSUBDIR}


if [ "${bBloreSumm}" = "true" ]; then
    BLORE_SUMMARY_JOBNAME="blore_${RANDSTRING}"
    for STUDY in ${STUDYNAMES[@]}; do
        JOBNAME="${BLORE_SUMMARY_JOBNAME}_${STUDY}"
        sed "s|_JOBNAME|${JOBNAME}|g;
             s|_B_LORE_|${BLORE}|g;
             s|_GEN_DIR|${THIS_LOCIDIR}/${STUDY}|g;
             s|_OUT_DIR|${BLORE_RESDIR}/${STUDY}|g;
             s|_SAM_FL_|${SAMPLEDIR}/${STUDY}/phenotypes.sample|g;
             s|_USELOCI|${USELOCI}|g;
            " ${MASTER_BSUBDIR}/blore_summary.bsub > ${JOBNAME}.bsub
        bsub < ${JOBNAME}.bsub
    done
fi

# Run B-LORE meta-analysis ========================================================
for NC in ${NCAUSAL}; do
    BLORE_META_JOBNAME="blore_meta_${NC}_${RANDSTRING}"
    sed "s|_JOBNAME|${BLORE_META_JOBNAME}|g;
         s|_B_LORE_|${BLORE}|g;
         s|_NCAUSAL|${NC}|g;
         s|_USELOCI|${USELOCI}|g;
         s|_OUT_DIR|${BLORE_METARS}|g;
         s|_STUDYN_|\"${STUDYNAMES[*]}\"|g;
         s|_BL_SUM_|${BLORE_RESDIR}|g;
        " ${MASTER_BSUBDIR}/blore_meta.bsub > ${BLORE_META_JOBNAME}.bsub
    if [ "${bBloreSumm}" = "true" ]; then
        bsub -w "done(${BLORE_SUMMARY_JOBNAME}*)" < ${BLORE_META_JOBNAME}.bsub
    else
        bsub < ${BLORE_META_JOBNAME}.bsub
    fi
done


cd ${THIS_JOBSUBDIR}
