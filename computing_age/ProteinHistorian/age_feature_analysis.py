#!/usr/bin/env python
"""
 age_feature_analysis.py - Copyright Tony Capra 2011

 Change log: 
 04/11/12 - updated to reflect paper revision
 01/18/12 - released


 This program supports the paper: 
 Capra JA, Williams AG, and Pollard KS. (2012) ProteinHistorian: Tools
 for the Comparative Analysis of Eukaryote Protein Origins. Submitted.

 Please cite the paper if you use this code.

 See usage and the README for more information and examples.

 This program was developed using Python 2.6.5, numpy 1.4.0, and
 matplotlib 1.0.0. Please contact the authors with any questions or
 bug reports.

 -----------------------------------------------------------------------------

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

"""

usage = """
USAGE:
age_feature_analysis.py [options] age_file feature_file

    -age_file, the ages for all proteins. Each line should have a one
     or more protein ids followed by a numeric age in the final
     column (tab-delimited).

    -feature_file, the protein features of interest. This file should
     have one protein id and the associated feature value on each
     line.  Any proteins in feature_file missing ages are ignored.

OPTIONS:
    -f [string]
     feature name. The name of the input feature data. Default=feature_file

    -h
     help. Print this message.

    -l [filename]
     age label file. A file mapping each protein with an age to string
     that is the "name" of the age. Default=None
        
    -n [string]
     name for output files. Base name given to all output
     files. Default=age_file--feature_file

    -o [path]
     output directory. Default=results/feature/output_name/

    -s 
     scale x-axis. Scale the x-axis of plots according to the ages;
     otherwise, each age bin is equally spaced when when printed. Default=False


"""

import os
import sys
import getopt
import datetime

import numpy
try:
    import scipy.stats as stats
except ValueError:
    import scipy.stats as stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pylab as plt


################################################################################

def parse_age_file(age_file, delimiter='\t'):
    """ Return prot2age and age2prot dictionaries and a list with the age
    of each protein in the file for use as a background. The age_file
    is tab-delimited with a variable number of ids for the same object
    and the age in the final column of each line. """

    prot2age = {}
    age2prot = {}
    background_ages = []

    for line in open(age_file):
	if line[0] == '#': continue

	t = line[:-1].split(delimiter)
        if len(t) < 2: continue

        prot_ids = t[:-1]
        age = None
        try:
            age = float(t[-1])
        except ValueError:
            age = t[-1]

        background_ages.append(age)

        for prot in prot_ids: 
            prot2age[prot] = age
            age2prot.setdefault(age, []).append(prot)

    return prot2age, age2prot, background_ages

def make_age2label_dict(age_label_file, age2prot):
    """ Return a dictionary mapping from ages found in age2prot to the
    corresponding entry for protein in age_label_file. It does not
    check that all proteins of a given age have the same label in the
    label file; it just takes the first label found. Usually these
    will be the names of the internal nodes corresponding to the depth
    ages, but they could, in theory, be anything."""

    age2label = {}
    # default label is the age itself
    for age in age2prot:
        age2label[age] = str(age)

    if age_label_file != None and os.path.exists(age_label_file):
        prot2label, label2prot, bg = parse_age_file(age_label_file)

        for age in age2prot:
            for prot in age2prot[age]:
                if age2label[age] != str(age): break

                if prot in prot2label:
                    age2label[age] = str(prot2label[prot]).replace('_', ' ')

    return age2label


################################################################################

def feature_by_age_boxplot(age2vals, age2label, outfn, title='', xlabel='Age',
                           ylabel='Value', scale_x = False, output_png=False, 
                           methods_str=''):
    fig = plt.figure()

    if len(age2vals.keys()) > 1:
        fig_width = fig.get_figwidth()
        fig_height = fig.get_figheight()
        fig.set_figwidth(fig_width * 1.7)
        fig.set_figheight(fig_height * 1.2)

    ax1 = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.26)
    
    ages = sorted(age2label.keys())

    box_data = []
    box_pos = []
    for i, age in enumerate(ages):
        if age in age2vals:
            box_data.append(age2vals[age])
            box_pos.append(age if scale_x else i)
        else:
            box_data.append([])
            box_pos.append(age if scale_x else i)

    bp = plt.boxplot(box_data, widths=.6, sym='', patch_artist=True, 
                     positions=box_pos)

    plt.setp(bp['boxes'], color='#99CCFF', edgecolor="black", lw=1)
    #plt.setp(bp['boxes'], color='darkkhaki', edgecolor="black", lw=1)
    plt.setp(bp['whiskers'], color='black', lw=1)
    plt.setp(bp['medians'], color='black', lw=1.5)
    plt.setp(bp['caps'], color='black', lw=1)


    for i, age in enumerate(ages):
        x_pos = age if scale_x else i
        if age in age2vals:
            ax1.plot(x_pos, numpy.average(age2vals[age]), 'x', color='red', 
                     markersize=6, markeredgewidth=1.5)


    if age2label: 
        labels = []
        for age in ages:
            if age < 0: 
                labels.append('')
            elif age in age2label:
                np = len(age2vals[age]) if age in age2vals else 0
                labels.append('%s (%d)' % (age2label[age], age))
                #labels.append('%s (n=%d)' % (age2label[age], np))
            else:
                labels.append('')

        xtickNames = plt.setp(ax1, xticklabels = labels)
        plt.setp(xtickNames, fontsize=10)
        plt.setp(xtickNames, rotation=45)
        plt.setp(xtickNames, horizontalalignment='right')


    ymin = ax1.viewLim.ymin
    ymax = ax1.viewLim.ymax
    y_range = ymax - ymin
    #ax1.set_ylim(-.05 * ymin, 1.05 * ymax)
    ax1.set_ylim(ymin - (.02 * y_range), ymax + (.02 * y_range))

    xmin = ax1.viewLim.xmin
    xmax = ax1.viewLim.xmax
    x_range = xmax - xmin
    pad = .02 * x_range
    ax1.set_xlim(xmin - pad, xmax + pad)


    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.title(title)

    plt.figtext(0.1, 0.01, '[DB: '+methods_str+']', size=10)

    plt.savefig(outfn)
    if output_png: plt.savefig(outfn.replace('.pdf', '.png'))
    
    return


################################################################################

def main():

    # parse command line options and update variables
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:l:n:o:sw")
    except getopt.GetoptError:
        sys.exit(usage)

    if len(args) != 2: 
        sys.exit(usage)

    AGE_FILE = args[0]
    FEATURE_FILE = args[1]

    FEATURE_NAME = FEATURE_FILE.split('/')[-1].split('.')[0]  # -f
    OUTPUT_DIR = 'results/feature/' # -o

    af_name = AGE_FILE.split('/')[-1]
    af_name = af_name.replace('.protein_list', '').replace('.txt', '')
    ff_name = FEATURE_FILE.split('/')[-1]
    ff_name = ff_name.replace('.protein_list', '').replace('.txt', '')
    OUTPUT_NAME = '%s--%s' % (af_name, ff_name) # -n 

    AGE_LABEL_FILE = None  # -l

    SCALE_X = False  # -s 

    WEB = False # -w, set to True if running on the web

    for opt, arg in opts:
        if opt == "-h":
            sys.exit(usage)

        elif opt == "-f":
            FEATURE_NAME = arg

        elif opt == "-l":
            AGE_LABEL_FILE = arg

        elif opt == "-n":
            OUTPUT_NAME = arg

        elif opt == "-o":
            OUTPUT_DIR = arg

        elif opt == "-s":
            SCALE_X = True

        elif opt == "-w":
            WEB = True


    if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)

    PREAMBLE = '# %s %s\n#\n' % (' '.join(sys.argv), datetime.datetime.now())
    PREAMBLE += '# OUTPUT_DIR = %s\n' % OUTPUT_DIR
    PREAMBLE += '# OUTPUT_NAME = %s\n' % OUTPUT_NAME
    PREAMBLE += '# FEATURE_NAME = %s\n' % FEATURE_NAME
    PREAMBLE += '# AGE_LABEL_FILE = %s\n' % AGE_LABEL_FILE
    PREAMBLE += '# SCALE_X = %s\n' % SCALE_X
    PREAMBLE += '#'

    # done parsing options

    # load data
    prot2age, age2prot, bg = parse_age_file(AGE_FILE)
    prot2feature, feature2prot, feature_bg = parse_age_file(FEATURE_FILE, 
                                                            delimiter=' ')
    # if no AGE_LABEL_FILE, labels are just ages
    if AGE_LABEL_FILE and not os.path.exists(AGE_LABEL_FILE): 
        AGE_LABEL_FILE = None
    age2label = make_age2label_dict(AGE_LABEL_FILE, age2prot)

    prots = []
    for prot in prot2feature:
        if prot in prot2age:
            prots.append(prot)

    PREAMBLE += '\n# %s: %d proteins.\n' % (FEATURE_FILE, len(prot2feature))
    PREAMBLE += '# %s: %d protein ids.\n' % (AGE_FILE, len(prot2age))
    PREAMBLE += '# Overlap: %d proteins.\n#' % (len(prots))

    # get list of ages and features and make age to length dictionary
    age2vals = {}
    obs_ages = []
    features = []
    for prot in prots:
        age2vals.setdefault(prot2age[prot], []).append(prot2feature[prot])
        obs_ages.append(prot2age[prot])
        features.append(prot2feature[prot])

    all_ages = sorted(age2prot.keys())

    # write stat file
    stat_filename = '%s%s_stats.txt' % (OUTPUT_DIR, OUTPUT_NAME)
    if not WEB: print 'writing %s...' % stat_filename
    stat_file = open(stat_filename, 'w')
    stat_file.write(PREAMBLE + '\n')

    corr, pval = stats.spearmanr(obs_ages,features)
    corr_str = 'Spearman rho=%.3g (p=%.3g), N=%d' % (corr, pval, len(obs_ages))
    stat_file.write('# ' + corr_str.replace('\n', '\n# ') + '\n')

    head_str = '# age\tlabel\tavg_%s\tmedian_%s\tnum_prot\n'
    stat_file.write( head_str % (FEATURE_NAME, FEATURE_NAME))
    for age in all_ages:
        if age in age2vals:
            data_str = '%s\t%s\t%.2f\t%.2f\t%d\n'
            stat_file.write( data_str % (age, age2label[age], 
                                         numpy.mean(age2vals[age]),
                                         numpy.median(age2vals[age]), 
                                         len(age2vals[age])))
        else:
            stat_file.write( '%s\t%s\tNone\tNone\t0\n' % (age, age2label[age]))
          
    stat_file.close()

    # make boxplot
    boxplot_filename = '%s%s_boxplot.pdf' % (OUTPUT_DIR, OUTPUT_NAME)
    if not WEB: print 'writing %s...' % boxplot_filename
    plot_title = 'Feature Distribution by Protein Age\n%s' % (corr_str)
    xlabel_str = 'Taxon of Origin'
    #if AGE_LABEL_FILE and '_age-time.' in AGE_LABEL_FILE: 
    if '_age-time.' in AGE_FILE: 
        xlabel_str = 'Taxon of Origin (age in millions of years)' 

    feature_by_age_boxplot(age2vals, age2label, boxplot_filename, 
                           title=plot_title, ylabel=FEATURE_NAME, 
                           xlabel=xlabel_str, output_png=WEB, 
                           scale_x=SCALE_X, methods_str=af_name)

    # End main()


if __name__ == "__main__": 
    main()
