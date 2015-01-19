#!/usr/bin/env python
"""
 age_enrichment_analysis.py - Copyright Tony Capra 2011

 Change log: 
 04/11/12 - updated to reflect paper revision
 01/18/12 - released


 This program supports the paper: 
 Capra JA, Williams AG, and Pollard KS. (2012) ProteinHistorian: Tools
 for the Comparative Analysis of Eukaryote Protein Origins. Submitted.

 Please cite the paper if you use this code.

 See usage and the README for more information and examples.

 Please contact the authors with any questions or bug reports.

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
age_enrichment_analysis.py [options] age_file protein_file1 [protein_file2]

    -age_file, the ages for all proteins. Each line should have a
     protein id followed by a numeric age (separated from the id by
     whitespace)

    -protein_file1, the proteins of interest. This file should
     list the identifiers of the proteins of interest one per line.  

    -protein_file2, an optional second list of proteins of interest to
     compare to the first. If not provided, the background of all
     proteins with ages in the species will be used.


OPTIONS:

    -a [string]
     handle mising poi. Sets how to handle an input protein of
     interest for which there is no age data. The options are: 'error'
     which stops the analysis; 'species_specific' which treats the
     protein as being specific to the species of interest; and
     'ignore' which ignores the protein. Default=error

    -h
     help. Print this message.

    -l [filename]
     age label file. A file mapping each protein with an age to string
     that is the "name" of the age. Default=None

    -m
     multiple hypothesis testing correction. Applies the Bonferroni
     correction to all tests.  Default=None
        
    -n [string]
     name for output files. Base name given to all output
     files. Default=age_file--feature_file

    -o [path]
     output directory. Default=results/enrichment/output_name/

    -p [string]
     name for proteins in protein_file1. A string giving a name for
     the first set of proteins of interest. Default=protein_file1

    -q [string]
     name for proteins in protein_file2. A string giving a name for
     the second set of proteins of interest. Default='Background' or
     protein_file2 if provided

    -s 
     scale x-axis. Scale the x-axis of plots according to the ages;
     otherwise, each age bin is equally spaced when when printed. This
     option only applies to the figure output; the input ages are
     still used in the analysis. Default=False

    -t [string]
     plot type. Either 'line' to output a line plot or 'bar' to output
     a bar plot. Default=bar

"""

import os
import sys
import getopt
import random
import datetime
import commands

import matplotlib
import matplotlib.pylab as plt
import seaborn as sns

import numpy
try:
    import scipy.stats as stats
except ValueError:
    import scipy.stats as stats



################################################################################

def parse_age_file(age_file):
    """ Return prot2age and age2prot dictionaries, a list with the age of
    each protein in the file for use as a background, and a dictionary
    mapping protein names to synonyms. The age_file is tab-delimited
    with a variable number of ids for the same object and the age in
    the final column of each line. """

    prot2age = {}
    age2prot = {}
    background_ages = []
    name2synonyms = {}

    for line in open(age_file):
	if line[0] == '#': continue

	t = line[:-1].split('\t')
        if len(t) < 2: continue

        prot_ids = t[:-1]
        age = None
        try:
            age = float(t[-1])
        except ValueError:
            age = t[-1]

        background_ages.append(age)

        for prot in prot_ids: 
            if prot == '': continue

            prot2age[prot] = age
            age2prot.setdefault(age, []).append(prot)
            name2synonyms[prot] = prot_ids

    return prot2age, age2prot, background_ages, name2synonyms

def parse_profile_file(phylo_profile_file):
    """ Return prot2profile, a dictionary mapping proteins to their phylo
    profile string, and species_str, a string with the species order for
    the profile. The phylo_profile_file should be generated by
    age_proteins.py. The first field gives the protein ids separated by a
    ~. """

    prot2profile = {}
    species_str = ''

    for line in open(phylo_profile_file):
	if line.startswith('# '): continue
        if line.startswith('#protein'): 
            species_str = line[:-1]

	t = line[:-1].split('\t')
        if len(t) < 2: continue

        prot_ids = t[0].split('~')

        for prot in prot_ids: 
            if prot == '': continue

            prot2profile[prot] = line[:-1]

    return prot2profile, species_str


def parse_poi(filename):
    """ Return a list of the proteins of interest. Protein ids can be on a
    single line separated by whitespace or on multiple lines or
    both. These should be checked against those in the input
    AGE_FILE."""

    poi = []
    for line in open(filename):
        if line[0] == '#': continue
        t = line[:-1].split()
        for prot_id in t:
            poi.append(prot_id)
        
    return poi

def parse_gaf_poi(filename, taxon_filter=None, evidence_filter=None):
    """ Return a list of lists in which the internal list contains
    possible IDs for each protein of interest. Potential IDs are taken
    from the DB Object ID, DB Object Symbol, DB Object Synonym
    columns. See http://www.geneontology.org/GO.format.gaf-2_0.shtml for
    the GAF 2.0 specification."""

    poi = []
    for line in open(filename):
        if line[0] == '!': continue
        t = line[:-1].split('\t')
        if len(t) < 16: continue

        taxon = t[12]
        evidence_code = t[6]

        if taxon_filter and taxon not in taxon_filter: continue
        if evidence_filter and evidence_code not in evidence_filter: continue

        poi.append([])
        poi[-1].append(t[2])  # name is more likely to match and filter dupls
        poi[-1].append(t[1])
        for id in t[10].split('|'): poi[-1].append(id)

    return poi

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
        prot2label, label2prot, bg, n2s = parse_age_file(age_label_file)

        for age in age2prot:
            for prot in age2prot[age]:
                if age2label[age] != str(age): break

                if prot in prot2label:
                    age2label[age] = str(prot2label[prot]).replace('_', ' ')
   
    return age2label


################################################################################

def age_plot(dataname2age2num, age2pval, age2label, outfn, plot_type='bar', 
             title='', output_png=False, scale_x=False, xlabel='Age', 
             methods_str=''):
    """  """
    ages = sorted(age2label.keys())

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    plt.subplots_adjust(bottom=0.27)

    width = 0.35

    plt.xlabel(xlabel)
    plt.ylabel('Fraction')

    # hack to make sure background is last if present
    data_names = []; bgs = []
    for name in sorted(dataname2age2num):
        if 'ackground' in name:
            bgs.append(name)
        else:
            data_names.append(name)
    data_names += bgs
      
    dataname2avgdepth = {}
    for dataname_idx, name in enumerate(data_names):
        # compute average depth, age, and total number of proteins
        N = 0.
        avg_depth = 0.
        avg_age = 0.
        for age in dataname2age2num[name]: 
            N += dataname2age2num[name][age]
            avg_depth += (age * dataname2age2num[name][age])

            try:
                avg_age += (float(age2label[age]) * dataname2age2num[name][age])
            except ValueError:
                avg_age += (age * dataname2age2num[name][age])

        avg_age /= N
        avg_depth /= N

        # plot data
        heights = []
        x_base = []
        for age_idx, age in enumerate(ages):
            if age in dataname2age2num[name] and \
                    dataname2age2num[name][age] != 0.0:
                heights.append(dataname2age2num[name][age] / N)
                x_base.append(age if scale_x else age_idx)
            else:
                heights.append(0.0)
                x_base.append(age if scale_x else age_idx)

        # line or bar plot
        if plot_type == 'line':
            line_color = '#990000' if dataname_idx % 2 == 0 else 'black'
            ax1.plot(x_base, heights, 'x-', label='%s (N=%d)' % (name, N), 
                     linewidth=2, mew=3, markersize=8, color=line_color)

        else:
            r = numpy.array(x_base, dtype='float_')
            bar_color = '#990000' if dataname_idx % 2 == 1 else 'black'
            r += .03 if dataname_idx % 2 == 1 else -.03

            ax1.bar(r - width * (1-dataname_idx), heights, width, 
                    label='%s (N=%d)' % (name, N), color=bar_color)



    # now set up axes limits and ticks
    ax1.set_xticks(x_base)
    x_range = x_base[-1] - x_base[0]
    pad = max(.03 * x_range, 1)
    ax1.set_xlim(x_base[0] - pad, x_base[-1] + pad)

    ymin = ax1.viewLim.ymin
    ymax = ax1.viewLim.ymax
    y_range = ymax - ymin

    ax1.set_ylim(-.05 * y_range, ymax + (.23 * y_range))

    leg = ax1.legend(prop={'size': 12}, loc='upper left', numpoints=1)
    leg.get_frame().set_fill(False)
    #leg.draw_frame(False)

    # write p-value *'s
    for age_idx, age in enumerate(ages):
        if age in age2pval:
            pval_str = ''
            if .01 < age2pval[age] < 0.05:
                pval_str = '*'
            elif .001 < age2pval[age] <= 0.01:
                pval_str = '**'
            elif age2pval[age] <= 0.001:
                pval_str = '***'

            ax1.text(age if scale_x else age_idx, ymin - (.025 * y_range), 
                     pval_str, horizontalalignment='center', 
                     verticalalignment='center', weight='bold')

    # set x_tick names
    labels = []
    for age in ages:
        if age in age2label:
            label = age2label[age]
            if label != str(age):
                labels.append('%s (%d)' % (label, age))
            else:
                labels.append('%d' % age)
        else:
            labels.append('')

    xtickNames = plt.setp(ax1, xticklabels = labels)
    plt.setp(xtickNames, fontsize=10)
    plt.setp(xtickNames, rotation=30)
    plt.setp(xtickNames, horizontalalignment='right')


    title_obj = plt.title(title, fontsize=12)

    plt.figtext(0.07, 0.01, '[DB: '+methods_str+']', size=9)

    plt.savefig(outfn)
    if output_png: plt.savefig(outfn.replace('.pdf', '.png'))
    
    return


################################################################################


def make_poi_age_list(poi_file, prot2age, handle_missing='error'):
    """ Return three lists. The names and ages (from prot2age) of the
    proteins in the set of interest and a list of proteins lacking
    ages. The approach to handling pois that are not found in the age
    dictionary is specified by handle_missing.  There are currently
    three approaches: 1) return an error and quit. 2) Ignore them. 3)
    Treat them as species specific proteins, i.e., having age 0. The
    poi_file can either have protein ids one per line, many per line
    separated by whitespace, or be in the GAF 2.0 format."""

    poi_names = []
    poi_ages = []
    poi_missing_ages = []

    f = open(poi_file)
    first_line = f.readline()
    f.close()

    if first_line.startswith('!gaf-version: 2.0'):
        poi = parse_gaf_poi(poi_file)
        for id_list in poi:
            selected_id = None
            for prot in id_list:
                if prot in prot2age: 
                    selected_id = prot
                    break

            if selected_id != None:
                if selected_id not in poi_names:
                    poi_names.append(selected_id)
                    poi_ages.append(prot2age[selected_id])

            else:
                if handle_missing == 'error':
                    sys.exit('ERROR: None of %s were found in the age file.' \
                                 % (', '.join(id_list)))
                elif handle_missing == 'species_specific':
                    poi_names.append(id_list[0])
                    poi_ages.append(0)
                    poi_missing_ages.append(id_list[0])
                elif handle_missing == 'ignore':
                    poi_missing_ages.append(id_list[0])
                    pass

    else: # parse simple ID list format
        poi = parse_poi(poi_file)
        for prot in poi:
            if prot in prot2age: 
                poi_names.append(prot)
                poi_ages.append(prot2age[prot])
            else:
                if handle_missing == 'error':
                    sys.exit('ERROR: %s was not found in the age file.' % (prot))
                elif handle_missing == 'species_specific':
                    poi_names.append(prot)
                    poi_ages.append(0)
                    poi_missing_ages.append(prot)
                elif handle_missing == 'ignore':
                    poi_missing_ages.append(prot)
                    pass

    return poi_names, poi_ages, poi_missing_ages


def fet_R(ct, alternative="t", filename='', get_ci=False):
    """ Use R to perform Fisher's exact test on the contingency table
    ct. Return the odds ratio and p-value.
    NOTE: This can be replaced with a call to stats.hypergeom.sf (1-cdf).
    sf(x,M,n,N) where x = # sub_type clustered, M = total # subs, n =
                          # sub_type, and N = # clustered subs
    """

    temp_file = filename
    if filename == '':
        temp_file = '%d.R' % random.randint(0, 10000000)

    f = open(temp_file, 'w') 
    f.write("ct <- matrix(c(%d, %d, %d, %d), nr = 2)\n" % (ct[0], ct[1], 
                                                           ct[2], ct[3]))
    f.write("fisher.test(ct, alternative = \"%s\")\n" % alternative)

    f.close()

    status, output = commands.getstatusoutput("Rscript %s" % temp_file)

    if filename == '':
        os.remove(temp_file)

    pval = None; odds_ratio = None; conf_int = None
    next_is_or = False
    next_is_ci = False
    for line in output.splitlines():
        if 'p-value' in line: 
            pval = float(line.split()[-1])
        elif line[:10] == 'odds ratio':
            next_is_or = True
        elif get_ci and line[:10] == '95 percent':
            next_is_ci = True
        elif next_is_or:
            next_is_or = False
            odds_ratio = float(line)
        elif get_ci and next_is_ci:
            next_is_ci = False
            t = line.split()
            conf_int = (float(t[0]), float(t[1]))
       
    return odds_ratio, pval, conf_int


def age_group_enrichment(poi1_ages, poi2_ages, 
                         poi2_not_bg=None, multi_correct=''):
    """ Perform Fisher's exact test on each age group.  The
    contingency table is membership in poi1 in one dimension and
    membership in the current age group in the other. Return
    dictionaries mapping ages to number in poi1, number in poi2,
    odds_ratios, and p-vals for difference. If multi_correct is not
    null, then apply the specified multiple test correction--currently
    only Bonferroni is implemented."""

    age2num_poi1 = {}; age2num_poi2 = {}  # one for each age
    age2pvals = {}; age2odds_ratio = {}
    ages = sorted(list(set(poi1_ages + poi2_ages)))
    for age in ages:
        # these are convenient for computing the cont. table
        num_poi1 = len(poi1_ages); num_poi2 = len(poi2_ages)
        num_poi1_age = poi1_ages.count(age)
        num_poi2_age = poi2_ages.count(age)

        c1 = num_poi1_age; c3 = num_poi1 - num_poi1_age
        # if poi2 is the background, then we have to subtract out poi1
        if poi2_not_bg == None:
            c2 = (num_poi2_age - num_poi1_age)
            c4 = (num_poi2 - num_poi2_age) - (num_poi1 - num_poi1_age)
        else: # we use the comparison set 
            c2 = num_poi2_age
            c4 = (num_poi2 - num_poi2_age)

        odds_ratio, pval, conf_int = fet_R((c1, c2, c3, c4))            
        #print c1, c2, c3, c4, odds_ratio, pval, conf_int

        age2num_poi1[age] = num_poi1_age
        age2num_poi2[age] = num_poi2_age
        age2odds_ratio[age] = odds_ratio
        age2pvals[age] = pval

    if multi_correct == 'bonferroni':
        num_test = len(age2pvals)
        for age in age2pvals:
            age2pvals[age] = min(age2pvals[age] * num_test, 1.0)

    return age2num_poi1, age2num_poi2, age2odds_ratio, age2pvals

################################################################################

def main():

    # Parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:hl:mn:o:p:q:st:wy:")
    except getopt.GetoptError:
        sys.exit(usage)

    if len(args) != 2 and len(args) != 3: 
        sys.exit(usage)

    AGE_FILE = args[0]
    POI1_FILE = args[1]
    POI2_FILE = None
    if len(args) == 3:
        POI2_FILE = args[2]

    POI1_NAME = '.'.join(POI1_FILE.split('/')[-1].split('.')[:-1]) # -p
    POI2_NAME = 'Background' # -q
    if POI2_FILE != None:
        POI2_NAME = '.'.join(POI2_FILE.split('/')[-1].split('.')[:-1]) 

    AGE_NAME = '.'.join(AGE_FILE.split('/')[-1].split('.')[:-1])
    OUTPUT_NAME = '%s--%s--%s' % (AGE_NAME, POI1_NAME, POI2_NAME) # -n
    OUTPUT_DIR = 'results/enrichment/%s/' % OUTPUT_NAME # -o

    PLOT_TYPE = 'bar'; PLOT_TYPES = ['line', 'bar'] # -t

    AGE_LABEL_FILE = None  # -l 
    PHYLO_PROFILE_FILE = None  # -y

    MULTI_CORRECT = False  # -m

    SCALE_X = False  # -s

    MISSING_POI = 'error' # -a 
    MISSING_POI_OPTIONS = ['error', 'species_specific', 'ignore']

    WEB = False  # -w, set if running on the web

    for opt, arg in opts:
        if opt == "-h":
            sys.exit(usage)

        elif opt == '-a':
            if arg in MISSING_POI_OPTIONS:
                MISSING_POI = arg
            else:
                err_str = "%s is not a valid missing poi option: %s" 
                sys.exit(err_str % (arg, ', '.join(MISSING_POI_OPTIONS)))

        elif opt == '-l':
            AGE_LABEL_FILE = arg

        elif opt == '-m':
            MULTI_CORRECT = 'bonferroni'

        elif opt == '-n':
            OUTPUT_NAME = arg

        elif opt == '-o':
            OUTPUT_DIR = arg

        elif opt == '-p':
            POI1_NAME = arg

        elif opt == '-q':
            POI2_NAME = arg

        elif opt == '-s':
            SCALE_X = True

        elif opt == '-t':
            if arg in PLOT_TYPES:
                PLOT_TYPE = arg
            else:
                err_str = "%s is not a valid plot type: %s"
                sys.exit(err_str % (arg, ', '.join(PLOT_TYPES)))

        elif opt == '-w':
            WEB = True

        elif opt == '-y':
            PHYLO_PROFILE_FILE = arg


    if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)

    PREAMBLE = '# %s %s\n#\n' % (' '.join(sys.argv), datetime.datetime.now())
    PREAMBLE += '# AGE_LABEL_FILE = %s\n' % AGE_LABEL_FILE
    PREAMBLE += '# PHYLO_PROFILE_FILE = %s\n' % PHYLO_PROFILE_FILE
    PREAMBLE += '# POI1_NAME = %s\n' % POI1_NAME
    PREAMBLE += '# POI2_NAME = %s\n' % POI2_NAME
    PREAMBLE += '# OUTPUT_NAME = %s\n' % OUTPUT_NAME
    PREAMBLE += '# OUTPUT_DIR = %s\n' % OUTPUT_DIR
    PREAMBLE += '#\n'
    PREAMBLE += '# MISSING_POI = %s\n' % MISSING_POI
    PREAMBLE += '# MULTI_CORRECT = %s\n' % MULTI_CORRECT
    PREAMBLE += '# PLOT_TYPE = %s\n' % PLOT_TYPE
    PREAMBLE += '# SCALE_X = %s\n' % SCALE_X
    PREAMBLE += '#\n'

    # done processing options

    # load data
    prot2age, age2prot, background_ages, name2syns = parse_age_file(AGE_FILE)
    all_ages = sorted(age2prot.keys())

    prot2profile = None
    if PHYLO_PROFILE_FILE: 
        prot2profile, species_order_str = parse_profile_file(PHYLO_PROFILE_FILE)

    # if no AGE_LABEL_FILE, labels are just ages
    if AGE_LABEL_FILE and not os.path.exists(AGE_LABEL_FILE): 
        AGE_LABEL_FILE = None
    age2label = make_age2label_dict(AGE_LABEL_FILE, age2prot)

    # significance analysis of proteins of interest
    poi1_names, poi1_ages, poi1_missing = make_poi_age_list(POI1_FILE, prot2age, 
                                                            handle_missing=MISSING_POI)

    poi2_names = None; poi2_missing = []  # don't print out the age file for poi2 if bg
    poi2_ages = background_ages   
    if POI2_FILE:
        poi2_names, poi2_ages, poi2_missing = make_poi_age_list(POI2_FILE, prot2age, 
                                                                handle_missing=MISSING_POI)

    # compute enrichment of poi in each age group
    age_info = age_group_enrichment(poi1_ages, poi2_ages, poi2_not_bg=POI2_FILE, 
                                    multi_correct=MULTI_CORRECT)
    age2num_poi1, age2num_poi2, age2odds_ratio, age2pvals = age_info

    # compute age summary stats
    avg_age1 = numpy.average(poi1_ages)
    median_age1 = numpy.median(poi1_ages)
    avg_age2 = numpy.average(poi2_ages)
    median_age2 = numpy.median(poi2_ages)

    # MWU Test for difference
    mwu, mwu_pval =  stats.mannwhitneyu(poi1_ages, poi2_ages)
    mwu_rho = mwu / (len(poi1_ages) * len(poi2_ages))


    ## Write output files

    # write age files
    outfn = '%s%s_poi1-ages.txt' % (OUTPUT_DIR, OUTPUT_NAME)
    outfile = open(outfn, 'w')
    if not WEB: print 'writing %s...' % outfn
    outfile.write(PREAMBLE)
    outfile.write('#protein\tage%s\n' % ('\tlabel' if AGE_LABEL_FILE else ''))
    for i, prot in enumerate(poi1_names):
        age = poi1_ages[i]
        label = ('\t' + age2label[age]) if AGE_LABEL_FILE else ''
        names = name2syns[prot] if prot in name2syns else prot
        outfile.write('%s\t%s%s\n' % ('~'.join(names), age, label))
    outfile.close()

    if PHYLO_PROFILE_FILE:
        profile_outfn = '%s%s_poi1-phylo-profiles.txt' % (OUTPUT_DIR, OUTPUT_NAME)
        profile_outfile = open(profile_outfn, 'w')
        if not WEB: print 'writing %s...' % profile_outfn
        profile_outfile.write(PREAMBLE)
        profile_outfile.write(species_order_str + '\n')
        for prot in poi1_names:
            if prot in prot2profile:
                profile_outfile.write('%s\n' % (prot2profile[prot]))
            else:
                profile_outfile.write('%s\tNot Found\n' % prot)
        profile_outfile.close()


    if poi2_names:
        outfn = '%s%s_poi2-ages.txt' % (OUTPUT_DIR, OUTPUT_NAME)
        outfile = open(outfn, 'w')
        if not WEB: print 'writing %s...' % outfn
        outfile.write(PREAMBLE)
        outfile.write('#protein\tage%s\n' % ('\tlabel' if AGE_LABEL_FILE else ''))
        for i, prot in enumerate(poi2_names):
            age = poi2_ages[i]
            label = ('\t' + age2label[age]) if AGE_LABEL_FILE else ''
            names = name2syns[prot] if prot in name2syns else prot
            outfile.write('%s\t%s%s\n' % ('~'.join(names), age, label))
        outfile.close()

        if PHYLO_PROFILE_FILE:
            profile_outfn = '%s%s_poi2-phylo-profiles.txt' % (OUTPUT_DIR, OUTPUT_NAME)
            profile_outfile = open(profile_outfn, 'w')
            if not WEB: print 'writing %s...' % profile_outfn
            profile_outfile.write(PREAMBLE)
            profile_outfile.write(species_order_str + '\n')
            for prot in poi1_names:
                if prot in prot2profile:
                    profile_outfile.write('%s\n' % (prot2profile[prot]))
                else:
                    profile_outfile.write('%s\tNot Found\n' % prot)
            profile_outfile.close()


    # write list of input proteins missing data
    missing_outfn = '%s%s_missing_pois.txt' % (OUTPUT_DIR, OUTPUT_NAME)
    if poi1_missing: 
        missing_outfile = open(missing_outfn, 'a')
        if MISSING_POI == 'ignore':
            print "WARNING: no ages found for %d proteins in %s:\n%s\n" % \
                (len(poi1_missing), POI1_NAME, ', '.join(poi1_missing))
            missing_outfile.write("WARNING: no ages found for %d proteins in %s:\n%s\n" % \
                (len(poi1_missing), POI1_NAME, ', '.join(poi1_missing)))

        elif MISSING_POI == 'species_specific':
            print "WARNING: considering %d proteins species specific due to lack of ages in %s:\n %s\n" % (len(poi1_missing), POI1_NAME, ', '.join(poi1_missing))
            missing_outfile.write("WARNING: considering %d proteins species specific due to lack of ages in %s:\n %s\n" % (len(poi1_missing), POI1_NAME, ', '.join(poi1_missing)))
        missing_outfile.close()

    if poi2_missing: 
        missing_outfile = open(missing_outfn, 'a')
        if MISSING_POI == 'ignore':
            print "WARNING: no ages found for %d proteins in %s:\n%s\n" % \
                (len(poi1_missing), POI1_NAME, ', '.join(poi1_missing))
            missing_outfile.write("WARNING: no ages found for %d proteins in %s:\n%s\n" % \
                (len(poi2_missing), POI2_NAME, ', '.join(poi2_missing)))

        elif MISSING_POI == 'species_specific':
            print "WARNING: considering %d proteins species specific due to lack of ages in %s:\n %s\n" % (len(poi1_missing), POI1_NAME, ', '.join(poi1_missing))
            missing_outfile.write("WARNING: considering %d proteins species specific due to lack of ages in %s:\n %s\n" % (len(poi2_missing), POI2_NAME, ', '.join(poi2_missing)))
        missing_outfile.close()


    # write results summary file
    outfn = '%s%s_summary.txt' % (OUTPUT_DIR, OUTPUT_NAME)
    outfile = open(outfn, 'w')
    if not WEB: print 'writing %s...' % outfn
    outfile.write(PREAMBLE)
    outfile.write('# Proteins in %s: %d\n' % (POI1_NAME, len(poi1_ages)))
    outfile.write('# Proteins in %s: %d\n#\n' % (POI2_NAME, len(poi2_ages)))

    outfile.write('# age\tlabel\t%s\t%s\tp-value\n' % (POI1_NAME, POI2_NAME))
    for age in all_ages:
        if age in age2num_poi1:  # if in poi1 also in poi2 by constr.
            format_str = '%s\t%s\t%d\t%d\t%.3f\n'
            outfile.write( format_str % (age, age2label[age], age2num_poi1[age], 
                                         age2num_poi2[age], age2pvals[age]) )
        else:
            outfile.write( '%s\t%s\t0\t0\tNone\n' % (age, age2label[age]) )


    mwu_str = 'Mann-Whitney U test: U = %.2g (p = %.3g)' % (mwu, mwu_pval)
    outfile.write('#\n# %s\n' % mwu_str)

    outfile.write('#\n# Average age of %s: %.1f\n' % (POI1_NAME, avg_age1))
    outfile.write('# Average age of %s: %.1f\n#\n' % (POI2_NAME, avg_age2))
    outfile.write('# Median age of %s: %.1f\n' % (POI1_NAME, median_age1))
    outfile.write('# Median age of %s: %.1f\n' % (POI2_NAME, median_age2))

    outfile.close()

    # Make plot(s)
    plot_name = '%s%s_%s.pdf' % (OUTPUT_DIR, OUTPUT_NAME, PLOT_TYPE)
    plot_title = 'Protein Age Distribution Comparison\n%s' % (mwu_str)

    xlabel_str = 'Taxon of Origin'
    if '-time.' in AGE_FILE: 
        xlabel_str = 'Taxon of Origin (age in millions of years)' 
        #xlabel_str = 'Age (millions of years ago)' 

    if not WEB: print 'writing %s...' % plot_name
    age_plot({POI1_NAME: age2num_poi1, POI2_NAME: age2num_poi2}, 
             age2pvals, age2label, plot_name, plot_type=PLOT_TYPE,
             title=plot_title, scale_x=SCALE_X, output_png=WEB, 
             xlabel=xlabel_str, methods_str=AGE_NAME)


    # End main()

if __name__ == "__main__": 
    main()
