import numpy as np
import collections
import os
import glob
import argparse

SNPINFO_FIELDS = ['chrm', 'rsid', 'bp_location', 'p']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()


def parse_args():

    parser = argparse.ArgumentParser(description='Filter SNPs for HWE and meta-presence after SNPTEST')

    parser.add_argument('-s', '--src',
                        type=str,
                        dest='srcdir',
                        metavar='DIR',
                        help='path of the directory containing loci definitions')

    parser.add_argument('-g', '--gen',
                        type=str,
                        dest='genofile',
                        metavar='FILE',
                        help='path of the genotypefile')

    parser.add_argument('-c', '--chr',
                        type=int,
                        dest='chrom',
                        metavar='CHROM',
                        help='chromosome number')

    parser.add_argument('-o', '--out',
                        type=str,
                        dest='outdir',
                        metavar='DIR',
                        help='path of the output directory')

    opts = parser.parse_args()
    return opts

opts = parse_args()

genotypefile = os.path.realpath(opts.genofile)
srcdir = os.path.realpath(opts.srcdir)
outdir = os.path.realpath(opts.outdir)

if not os.path.exists(outdir):
    os.makedirs(outdir)

nloci = len(glob.glob1(srcdir, "*.selected"))

loci = [[] for i in range(nloci)]
for i in range(nloci):
    locus = i + 1
    infile = os.path.join(srcdir, "Locus.{:03d}.selected".format(locus))
    with open(infile, 'r') as rfile:
        next(rfile)
        next(rfile)
        for line in rfile:
            linesplit = line.split()
            chrm = int(linesplit[0])
            rsid = linesplit[1]
            bp = int(linesplit[2])
            pval = float(linesplit[3])
            snp = SnpInfo(chrm = chrm, rsid = rsid, bp_location = bp, p = pval)
            loci[i].append(snp)

chrmloci = [x for x in loci if x[0].chrm == opts.chrom]
mlines = collections.defaultdict(lambda:0)
rsidlist = list()

if len(chrmloci) > 0:
    # Read each chromosome and store in memory only if there is one or more loci
    with open(genotypefile, 'r') as rfile:
        for line in rfile:
            rsid = line.split()[1]
            rsidlist.append(rsid)
            mlines[rsid] = line

for locus in chrmloci:
    indx = loci.index(locus) + 1
    outfile = os.path.join(outdir, "Locus.{:03d}.gen".format(indx))
    with open(outfile, 'w') as wfile:
        for snp in locus:
            if snp.rsid in rsidlist:
                wfile.write(mlines[snp.rsid])
