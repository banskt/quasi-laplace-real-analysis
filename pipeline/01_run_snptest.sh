#!/bin/bash

source PATHS

RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
THIS_SIMDIR="${BASEDIR}/metaanalysis"
SNPTEST_JOBSUBDIR="${JOBSUBDIR}/metaanalysis"

for STUDY in ${STUDYNAMES[@]}; do
    SAMPLEDIR="${THIS_SIMDIR}/samples/${STUDY}"
    if [ ! -d ${SAMPLEDIR} ]; then mkdir -p ${SAMPLEDIR}; fi
    THIS_SAMPLE=${SAMPLE_FMT//\[STUDY\]/${STUDY}}
    cp ${THIS_SAMPLE} ${SAMPLEDIR}/phenotypes.sample
done

if [ ! -d ${SNPTEST_JOBSUBDIR} ]; then mkdir -p ${SNPTEST_JOBSUBDIR}; fi
cd ${SNPTEST_JOBSUBDIR}

SNPTEST_JOBNAME="snptest_${RANDSTRING}"
for STUDY in ${STUDYNAMES[@]}; do
    _GENO_STUDY_FMT=${GENO_FMT//\[STUDY\]/${STUDY}}
    for CHROM in {1..22}; do
        JOBNAME="${SNPTEST_JOBNAME}_${STUDY}_${CHROM}"
        GENOTYPEFILE=${_GENO_STUDY_FMT//\[CHRM\]/${CHROM}}
        sed "s|_JOBNAME|${JOBNAME}|g;
             s|_SIMDIR_|${THIS_SIMDIR}|g;
             s|_GSTUDY_|${STUDY}|g;
             s|_SNPTEST|${SNPTEST}|g;
             s|_USEGENO|${GENOTYPEFILE}|g;
             s|_CHROM__|${CHROM}|g;
            " ${MASTER_BSUBDIR}/snptest.bsub > ${JOBNAME}.bsub
        bsub < ${JOBNAME}.bsub
    done
done

META_JOBNAME="meta_${RANDSTRING}"
sed "s|_JOBNAME|${META_JOBNAME}|g;
     s|_SIMDIR_|${THIS_SIMDIR}|g;
     s|_STUDYN_|\"${STUDYNAMES[*]}\"|g;
     s|_SAMPLES|\"${STUDYSAMPLES[*]}\"|g;
     s|_GENINF_|${GENINF}|g;
     s|_FILSNP_|${FILTERSNPS}|g;
     s|__META__|${META}|g;
    " ${MASTER_BSUBDIR}/meta.bsub > ${META_JOBNAME}.bsub
bsub -w "done(${SNPTEST_JOBNAME}*)" < ${META_JOBNAME}.bsub

cd ${CURDIR}
