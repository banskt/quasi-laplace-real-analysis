import numpy as np
import scipy.stats
import collections
import os
import argparse
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
onemb = 1000000

def parse_args():

    parser = argparse.ArgumentParser(description='Select the loci from meta-analysis')

    parser.add_argument('-d', '--resdir',
                        type=str,
                        dest='resdir',
                        metavar='DIR',
                        help='path of the directory containing results of SNPTEST/META')

    parser.add_argument('-p', 
                        type=float,
                        dest='pcut',
                        metavar='FLOAT',
                        help='Minimum p-value used for cut-off')

    parser.add_argument('-c', '--criteria',
                        type=str,
                        dest='criteria',
                        metavar='STR',
                        help='the criteria for loci creation -- how many SNPs to include -- any number or \'all\'')


    parser.add_argument('-o', '--out',
                        type = str,
                        dest = 'outdir',
                        metavar = 'DIR',
                        help = 'path of the output directory')

    parser.add_argument('-k', '--known',
                        type = str,
                        dest = 'known',
                        metavar = 'FILE',
                        help = 'path of the known loci')

    opts = parser.parse_args()
    return opts


def second_smallest(numbers):
    count = 0
    m1 = m2 = float('inf')
    for x in numbers:
        count += 1
        if x < m2:
            if x <= m1:
                m1, m2 = x, m1            
            else:
                m2 = x
    return m2 if count >= 2 else None


SNPINFO_FIELDS = ['chrm', 'rsid', 'bp_location', 'p']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()

LOCUSINFO_FIELDS = ['chrm', 'rsid', 'bp_location', 'gene']
class LocusInfo(collections.namedtuple('_LocusInfo', LOCUSINFO_FIELDS)):
    __slots__ = ()

opts = parse_args()
pcut = opts.pcut
bprange = int(0.4 * onemb)
nchrm = 22
criteria = opts.criteria
known_loci_file = os.path.realpath(opts.known)
resdir = os.path.realpath(opts.resdir)
outdir = os.path.realpath(opts.outdir)
if not os.path.exists(outdir):
    os.makedirs(outdir)
plotfile = os.path.join(outdir, 'selected_loci.png')

# Read all the SNP information from meta-analysis
print ('Reading meta-analysis results')
snpinfo = list()
snpinfo_chrm = [[] for i in range(nchrm)]
for i in range(nchrm):
    chrm = i + 1
    infile = os.path.join(resdir, 'chr{:d}_qc_meta_fixed_effect.out'.format(chrm))
    with open(infile, 'r') as rfile:
        next(rfile)
        for line in rfile:
            linesplit = line.split()
            rsid = linesplit[1]
            bp = int(linesplit[2])
            pval = float(linesplit[5])
            snp = SnpInfo(chrm = chrm, rsid = rsid, bp_location = bp, p = pval)
            snpinfo.append(snp)
            snpinfo_chrm[i].append(snp)

# Read the position of the known loci
print('Reading known loci positions')
known_loci = list()
with open(known_loci_file, 'r') as rfile:
    next(rfile)
    for line in rfile:
        linesplit = line.split()
        chrm = int(linesplit[0])
        rsid = linesplit[1]
        bppos = int(linesplit[2])
        gene = linesplit[4]
        locus = LocusInfo(chrm = chrm, rsid = rsid, bp_location = bppos, gene = gene)
        known_loci.append(locus)

#Select the lead SNPS of the loci
print ('Selecting lead SNPs')
pvals = np.array([snp.p for snp in snpinfo])
pord = np.argsort(pvals)
leadsnps = list()
i = 0
#while len(leadsnps) < nloci:
while pvals[pord[i]] <= pcut:
    index = pord[i]
    snp = snpinfo[index]
    chrm = snp.chrm
    # Check previous SNPs from the same chromosome
    prev_leads = [x for x in leadsnps if x[0].chrm == chrm]
    #iprev = [i for i in range(len(leadsnps)) if leadsnps[i][0].chrm == chrm]
    iprev = [leadsnps.index(x) for x in prev_leads]
    if len(prev_leads) == 0:
        # no SNPs from the chromosome. Append.
        leadsnps.append([snp])
    else:
        # calculate distance from previous lead SNPs in the same chromosome
        dists = [min([abs(snp.bp_location - x.bp_location) for x in sublist]) for sublist in prev_leads]
        #print(dists)
        if min(dists) > bprange:
            # Well separated from previous lead snps in the chromosome. Append.
            leadsnps.append([snp])
        else:
            k = dists.index(min(dists))
            indx = iprev[k]
            other_leads = [x for x in prev_leads if not x == leadsnps[indx]]
            iother = [leadsnps.index(x) for x in other_leads]
            leadsnps[indx].append(snp)
            if len(other_leads) > 0:
                dists = [min([abs(snp.bp_location - x.bp_location) for x in sublist]) for sublist in other_leads]
                if min(dists) < bprange:
                    # Merge the two SNPs
                    k = dists.index(min(dists))
                    nindx = iother[k]
                    leadsnps[indx] += leadsnps[nindx]
                    leadsnps.pop(nindx)
    i += 1

nloci = len(leadsnps)

# Create the loci for all lead SNPs
# In the first loop, we do not filter SNPs for p-values, we select all of them.
print ('Creating the loci')
locilist = [[] for i in range(nloci)]
info = ['None' for i in range(nloci)]
causality = [0 for i in range(nloci)]
mindist_lead = [None for i in range(nloci)]
unused = [[] for i in range(nloci)]
diff = int(bprange / 2)
realgenecounter = 0
reallocuscounter = 0
for i in range(nchrm):
    chrm = i + 1
    leads = [x for x in leadsnps if x[0].chrm == chrm]
    nlead = len(leads)
    ileads = [leadsnps.index(x) for x in leads]
    bp = [[snp.bp_location for snp in lead] for lead in leads]
    leadrange = [range(min(x) - diff, max(x) + diff ) for x in bp]
    for snp in snpinfo_chrm[i]:
        indx = [ileads[i] for i in range(nlead) if snp.bp_location in leadrange[i]]
        if len(indx) > 1:
            print("Overlapping ranges")
        elif len(indx) == 1:
            locilist[indx[0]].append(snp)
        elif len(indx) == 0:
            unused[i].append(snp)

    c4dloci = [x for x in known_loci if x.chrm == chrm]
    for j in range(nlead):
        if len(c4dloci) > 0:
            a = [x for x in c4dloci if x.bp_location in leadrange[j]]
            ngene = len(a)
            if ngene > 0:
                rsid = ''
                gene = ''
                for k in range(ngene):
                    rsid += '%s ' % a[k].rsid
                    gene += '%s ' % a[k].gene
                causality[ileads[j]] = 1
                mindist_lead[ileads[j]] = 0
                info[ileads[j]] = 'Causality: 1 // Lead SNP(s) %s, gene(s) %s' % (rsid, gene)
                realgenecounter += len(set([x.gene for x in a]))
                reallocuscounter += 1
            else:
                dists = [min(abs(x.bp_location - leadrange[j][0]), abs(x.bp_location - leadrange[j][1])) for x in c4dloci]
                causality[ileads[j]] = 0
                mindist_lead[ileads[j]] = min(dists) / onemb
                info[ileads[j]] = 'Causality: 0 // Distance from nearest lead SNP is %g Mb' % mindist_lead[ileads[j]]
        else:
            causality[ileads[j]] = 0
            mindist_lead[ileads[j]] = None
            info[ileads[j]] = 'Causality: 0 // No locus in this chromosome.'

locilist_filter = list()
info_filter = list()
causality_filter = list()
mindist_filter = list()
# Filter out loci with single abnormally high p-value SNP
print ('Filtering the loci')
for i, locus in enumerate(locilist):
    print ("Locus %i with %i SNPs" % (i + 1, len(locus)))
    pvals = [x.p for x in locus]
    pval2 = second_smallest(pvals)
    if pval2 > 0.005:
        # remove this loci
        if causality[i] == 1:
            print ("Removing causal locus %i from chr%i, %g Mb" % (i+1, locus[0].chrm, locus[0].bp_location / onemb))
            print (info[i], pval2)
        else:
            print ("Removing locus %i from chr%i, %g Mb" % (i+1, locus[0].chrm, locus[0].bp_location / onemb))
            print (info[i], pval2)
        rejected = [x for x in locus]
        chrm = locus[0].chrm
        unused[chrm - 1] += rejected
    else:
        # keep the locus
        locilist_filter.append(locus)
        causality_filter.append(causality[i])
        mindist_filter.append(mindist_lead[i])
        info_filter.append(info[i])

locilist = locilist_filter
info = info_filter
mindist_lead = mindist_filter
causality = causality_filter

#newloci = locilist
# If required, apply a filter of p-values to select SNPs per locus
print ('Selecting %s SNPs' % criteria)
if criteria == 'all':
    newloci = locilist
    #newloci = locilist[:25]
    #rejectedloci = locilist[25:]
    #for i, locus in enumerate(rejectedloci):
    #    rejected = [x for x in locus]
    #    chrm = locus[0].chrm
    #    unused[chrm - 1] += rejected
else:
    maxtarget = int(criteria)
    oldloci = locilist
    newloci = list()
    for i, locus in enumerate(oldloci):
        #nsnps = len(locus)
        target = min(len(locus), maxtarget)
        #nlead = len(leadsnps[i])
        #if nlead < 50:
        #    target = 100
        #else:
        #    target = nlead * 2
        plist = np.array([x.p for x in locus])
        kpval = list(np.argsort(plist)[:target])
        bplist = [locus[i].bp_location for i in kpval]
        kbpos = sorted(range(target), key=lambda k: bplist[k])
        newlocus = [locus[kpval[i]] for i in kbpos]
        newloci.append(newlocus)
        rejected = [x for x in locus if x not in newlocus]
        chrm = locus[0].chrm
        unused[chrm - 1] += rejected



# Lets see how the loci looks.
print ('Plotting the selections')
kelly_colors_hex = [
    #'#FFB300', # Vivid Yellow
    '#803E75', # Strong Purple
    '#FF6800', # Vivid Orange
    '#A6BDD7', # Very Light Blue
    '#C10020', # Vivid Red
    '#CEA262', # Grayish Yellow
    '#817066', # Medium Gray

    # The following don't work well for people with defective color vision
    '#007D34', # Vivid Green
    '#F6768E', # Strong Purplish Pink
    '#00538A', # Strong Blue
    '#FF7A5C', # Strong Yellowish Pink
    '#53377A', # Strong Violet
    '#FF8E00', # Vivid Orange Yellow
    '#B32851', # Strong Purplish Red
    '#F4C800', # Vivid Greenish Yellow
    '#7F180D', # Strong Reddish Brown
    '#93AA00', # Vivid Yellowish Green
    '#593315', # Deep Yellowish Brown
    '#F13A13', # Vivid Reddish Orange
    '#232C16', # Dark Olive Green
    ]
colors = kelly_colors_hex #['red' for i in range(200)]
figsize = (40, 15)
fig = plt.figure(figsize=figsize)
nrow = 4
ncol = 6
for i in range(nchrm):
    ax = fig.add_subplot(nrow, ncol, i+1)
    chrm = i + 1
    loci = [x for x in newloci if x[0].chrm == chrm]
    bpmax = snpinfo_chrm[chrm - 1][-1].bp_location
    bpmax = round(bpmax / 1000000, 0)
    bp = [x.bp_location / 1000000 for x in unused[i]]
    p = [-np.log10(x.p) for x in unused[i]]
    ax.scatter(bp, p, color='gainsboro', s = 2, alpha=0.1)
    for i, locus in enumerate(loci):
        #print(i)
        bp = [x.bp_location / 1000000 for x in locus]
        p = [min(-np.log10(x.p), 9.8) for x in locus]
        ax.scatter(bp, p, color = colors[i % len(colors)], s = 10, alpha = 0.6)
    ax.set_xlim(0, bpmax)
    ax.set_ylim(-1, 10)
    ax.set_title("Chr%s" % chrm)
plt.tight_layout()
#plt.show()
plt.savefig(plotfile)

# Finally, save the loci
print ('Creating output')
for i, locus in enumerate(newloci):
    outfile = os.path.join(outdir, 'Locus.{:03d}.selected'.format(i + 1))
    with open(outfile, 'w') as wfile:
        wfile.write("%s\n" % info[i])
        wfile.write("Chr rsid bp p-value\n")
        for snp in locus:
            wfile.write("%i %s %i %g\n" % (snp.chrm, snp.rsid, snp.bp_location, snp.p))

print ('%i locus created with %g Mb flanking. %i causal loci. Before filter: %i causal loci containing %i genes.' % (len(newloci), diff / onemb, sum(causality), reallocuscounter, realgenecounter))
