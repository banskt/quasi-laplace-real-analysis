import os
import numpy as np
import scipy.stats
import argparse
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

def parse_args():

    parser = argparse.ArgumentParser(description='Calculate genomic inflation factor in the simulations')

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

    parser.add_argument('-o', '--outfile',
                        default='genomic_inflation_factor.txt',
                        type=str,
                        dest='outfilename',
                        metavar='FILE',
                        help='name of output file for genomic inflation factors')

    opts = parser.parse_args()
    return opts


opts = parse_args()
snptestdir = os.path.realpath(opts.snptestdir)
studies = opts.studies
outfilename = os.path.realpath(opts.outfilename)

tol = 1e-5

outfile = open(outfilename, 'w')
for study in studies:
    pvalsfile = os.path.join(snptestdir, study, 'pvals.dat')
    pvals = np.loadtxt(pvalsfile)
    ordered_pvals = pvals[np.argsort(pvals)]
    exp_pvals = np.random.uniform(0, 1, pvals.shape[0])
    ordered_exp_pvals = exp_pvals[np.argsort(exp_pvals)]
    logpval = - np.log10(ordered_pvals)
    logpval_exp = - np.log10(ordered_exp_pvals)
    x = np.arange(0, 10)
    y = x

    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.scatter(logpval_exp, logpval)
    ax1.plot(x, y)
    ax1.set_ylim([0,10])
    ax1.set_xlim([0,10])
    figfile = os.path.join(snptestdir, study, 'qqplot.pdf')
    plt.savefig(figfile)

    sel_pvals = pvals[np.where(pvals > tol)]
    print (sel_pvals.shape)
    chisq = scipy.stats.chi2.ppf(1 - sel_pvals, 1)
    geninf = np.median(chisq) / scipy.stats.chi2.ppf(0.5, 1)

    outfile.write("%s %s\n" % (study, geninf))
outfile.close()
