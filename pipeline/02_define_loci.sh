#!/bin/bash

source PATHS
source CONFIG

RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
METARESDIR="${BASEDIR}/metaanalysis/snptest/meta"
THIS_JOBSUBDIR="${JOBSUBDIR}/${CONFIGDIR}/selectloci/"
THIS_CONFIGDIR="${BASEDIR}/${CONFIGDIR}/selectloci"

if [ ! -d ${THIS_JOBSUBDIR} ]; then mkdir -p ${THIS_JOBSUBDIR}; fi
cd ${THIS_JOBSUBDIR}

LOCI_JOBNAME="define_loci_${RANDSTRING}"
LOCIDEFDIR="${THIS_CONFIGDIR}/locidef"
sed "s|_JOBNAME|${LOCI_JOBNAME}|g;
     s|_OUTDIR_|${LOCIDEFDIR}|g;
     s|_METDIR_|${METARESDIR}|g;
     s|_DEFLOC_|${DEFINELOCI}|g;
     s|_PVL_MN_|${PVAL_CUTOFF_LOCI_SELECTION}|g;
     s|_MX_SNP_|${MAX_SNPS_EACH_LOCUS}|g;
     s|_KN_LOC_|${KNOWNLOCI}|g;
    " ${MASTER_BSUBDIR}/defineloci.bsub > ${LOCI_JOBNAME}.bsub

bsub < ${LOCI_JOBNAME}.bsub

CREATE_JOBNAME="create_loci_${RANDSTRING}"
for STUDY in ${STUDYNAMES[@]}; do
    _GENO_STUDY_FMT=${GENO_FMT//\[STUDY\]/${STUDY}}
    OUTDIR="${THIS_CONFIGDIR}/${STUDY}"
    for CHROM in {1..22}; do
        THIS_JOBNAME="${CREATE_JOBNAME}_${STUDY}_${CHROM}"
        GENOTYPEFILE=${_GENO_STUDY_FMT//\[CHRM\]/${CHROM}}
        sed "s|_JOBNAME|${THIS_JOBNAME}|g;
             s|_CRTLOC_|${CREATELOCI}|g;
             s|_LOCDEF_|${LOCIDEFDIR}|g;
             s|_USEGENO|${GENOTYPEFILE}|g;
             s|_CHR_NM_|${CHROM}|g;
             s|_OUTDIR_|${OUTDIR}|g;
            " ${MASTER_BSUBDIR}/createloci.bsub > ${THIS_JOBNAME}.bsub
        bsub -w "done(${LOCI_JOBNAME}*)" < ${THIS_JOBNAME}.bsub
    done
done

cd ${CURDIR}
