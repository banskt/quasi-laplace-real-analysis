import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import collections
import argparse
import os

def parse_args():

    parser = argparse.ArgumentParser(description='Filter SNPs for HWE and meta-presence after SNPTEST')

    parser.add_argument('-d', '--snptestdir',
                        type=str,
                        dest='snptestdir',
                        metavar='DIR',
                        help='path of the directory containing results of SNPTEST')

    parser.add_argument('-s', '--studies',
                        nargs='*',
                        type=str,
                        dest='studies',
                        metavar='STR',
                        help='list of study names')

    opts = parser.parse_args()
    return opts

opts = parse_args()
snptestdir = os.path.realpath(opts.snptestdir)
studies = opts.studies

nstudy = len(studies)
phwe_tol = 1e-4
nstudy_tol = 3
for i in range(22):
    chrm = i + 1
    print ("Selecting SNPs of chromosome %i." % chrm)
    presentin = collections.defaultdict(lambda:0)
    phwe = [collections.defaultdict(lambda:0) for i in range(nstudy)]
    greaterphwe = collections.defaultdict(lambda:0)
    
    for j, study in enumerate(studies):
        hwefile = os.path.join(snptestdir, "{:s}/chr{:d}.hwe".format(study, chrm))
        with open(hwefile, 'r') as rfile:
            for line in rfile:
                if not line.startswith('#') and not line.startswith('alternate_ids'):
                    linesplit = line.split()
                    rsid = linesplit[1]
                    phwe[j][rsid] = float(linesplit[20])
    
    for study in studies:
        snptestfile = os.path.join(snptestdir, "{:s}/chr{:d}.out".format(study, chrm))
        with open(snptestfile, 'r') as rfile:
            for line in rfile:
                if not line.startswith('#') and not line.startswith('alternate_ids'):
                    linesplit = line.split()
                    rsid = linesplit[1]
                    presentin[rsid] += 1
                    
    for rsid, val in presentin.items():
        a = [phwe[k][rsid] > phwe_tol for k in range(nstudy)]
        greaterphwe[rsid] = sum(a)

    for k, study in enumerate(studies):
        nsnps = 0
        nqc = 0
        nmeta = 0
        snptestfile = os.path.join(snptestdir, "{:s}/chr{:d}.out".format(study, chrm))
        outfile     = os.path.join(snptestdir, "{:s}/chr{:d}_qc.out".format(study, chrm))
        
        with open(outfile, 'w') as wfile:
            with open(snptestfile, 'r') as rfile:
                for line in rfile:
                    if not line.startswith('#') and not line.startswith('alternate_ids'):
                        linesplit = line.split()
                        rsid = linesplit[1]
                        nsnps += 1
                        if presentin[rsid] >= nstudy_tol:
                            nmeta += 1
                            if greaterphwe[rsid] >= nstudy_tol and phwe[k][rsid] > phwe_tol:
                                nqc += 1
                                wfile.write(line)
                    else:
                        wfile.write(line)
                        
        print ("Chr%i %s: Total %i SNPs, after QC %i SNPs. %i SNPs failed meta-presence, %i SNPs failed HWE"
               % (chrm, study, nsnps, nqc, nsnps - nmeta, nmeta - nqc))

