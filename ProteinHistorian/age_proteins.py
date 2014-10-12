#!/usr/bin/env python
"""
 age_proteins.py - Copyright Tony Capra 2011

 Change log: 
 04/11/12 - updated to include wagner parsimony
 01/18/12 - released

 This program supports the paper: 
 Capra JA, Williams AG, and Pollard KS. (2012) ProteinHistorian: Tools
 for the Comparative Analysis of Eukaryote Protein Origins. Submitted.

 Please cite the paper if you use this code.

 See usage and the README for more information and examples. This
 program was developed using Python 2.6.5. It also requires DendroPy
 3.7.0 and the Count package
 (http://www.iro.umontreal.ca/~csuros/gene_content/count.html). Please
 contact the authors with any questions or bug reports.

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
age_proteins.py [options] family_file species_tree

    -family_file, the partition of all proteins in all species into
     evolutionarily related families. Each protein from each species
     should appear once and only once. In this file, each protein is
     identified by its species, a divider token, and its
     identifier. For example: HUMAN|ENSG00000143858 MOUSE|MGI:99666
     RAT|3804

    -species_tree, The tree must be in the newick format and each leaf
     node must have a unique label.  If internal nodes of the tree are
     labeled, then these labels will be used as possible ages. If not,
     then the leaf nodes under a given internal node will be
     concatenated to provide a name for the node.

OPTIONS:

    -a [string]
     ancestral reconstruction algorithm. ProteinHistorian provides two
     algorithms for reconstructing the ancestral history of each
     protein family: Dollo parismony ('dollo') and asymmetric Wagner
     parsimony ('wagner'). Wagner parsimony analysis requires the
     Count program (Csuros 2010). A gain penalty can be set (-g) to
     weight gains and losses in the Wagner parsimony analysis. We
     recommend using Wagner parsimony.  Default='dollo'

    -f [comma separated string]
     output formats. Sets the types of age files output. The options
     are: 'depth' in which each protein's age is given in terms of
     depth in the species tree; 'label' outputs the label of the node
     of origin as the age; 'time' uses divergence estimates from
     TimeTree to assign ages in millions of years ago; and 'profile'
     prints a phylogenetic profile for each protein.
     Default='depth,label,time,profile'

    -g [positive real number]
     gain penalty. The relative penalty for gains compared to losses
     in Wagner parsimony. g < 1 favors gains, while g > 1 favors
     losses.  Default=1

    -h
     help. Print this message.

    -k [character]
     species name delimiter. This character is taken as the delimiter
     between the spcies name and protein id in the family file,
     e.g. HUMAN|ENSG00000143858. Default=|

    -n [string]
     output base name. Base name given to all output
     files. Default=family_file prefix

    -o [path]
     output directory. Default=results/ages/output_name/

    -s [comma separated string]
     species. Sets the species for which to output age files. If no
     species are given, then ages are output for all species in the
     input tree. Default=''

    -t [path]
     node age file. Path to a file that maps node names in the species
     tree to estimated divergence times. Default=None

"""

import os
import sys
import getopt
import random
import commands
import datetime

import dendropy


################################################################################


def parse_families(filename, valid_names = None, species_token='|'):
    """ Parse the families defined in filename. Return a dictionary
    mapping sequence names to families, families to sequence names,
    and species names to sequence names. If valid_names is provided,
    only include seqs listed. In the file each family is represented
    by a single line, and its members are separated by whitespace.
    The name of each protein must have the following form:
    species_id|protein_id, where | is species_token."""
    prot2family = {}
    family2prots = {}
    species2prots = {}

    n = 0
    for line in open(filename):
        if line[0] == '#': continue
        t = line[:-1].split()

        n += 1
        family_name = 'family%d' % n

        family2prots[family_name] = []
        for prot in t:
            if valid_names == None or prot in valid_names:
                species = prot.split(species_token)[0]
                species2prots.setdefault(species, []).append(prot)

                prot2family[prot] = family_name
                family2prots[family_name].append(prot)
              
    return prot2family, family2prots, species2prots


def write_attribute_file(outfn, attr_dict, preamble=None):
    """ Write each key, value pair in attr_dict to outfn. """

    outfile = open(outfn, 'w')
    if preamble != None: outfile.write(preamble)

    for key in sorted(attr_dict):
        outfile.write('%s\t%s\n' % (key, attr_dict[key]))

    return outfile.close()

def write_attribute_file_no_species(outfn, attr_dict, preamble=None, 
                                    combine_ids=False, additional_id_dict=None, 
                                    species_token='|'):
    """ Write each key, value pair in attr_dict to outfn. But strip the
    species ID and token from each key. If combine_ids, print all ids
    using combine_ids as the separator. This is used for the phylo
    profile. If extra_id_dict exists, check each protein name for
    additional names to print."""

    outfile = open(outfn, 'w')
    if preamble != None: outfile.write(preamble)

    for key in sorted(attr_dict):
        # this turns HUMAN|DB1:ID1|DB2:ID2 into ['ID1', 'ID2']
        names = [':'.join(x.split(':')[1:]) for x in key.split(species_token)[1:]]

        print_names = names[:]        
        if additional_id_dict:
            for name in names:
                if name in additional_id_dict:
                    print_names += additional_id_dict[name]

        if combine_ids:
            outfile.write('%s\t%s\n' % (str(combine_ids).join(print_names), 
                                        attr_dict[key]))
        else:
            outfile.write('%s\t%s\n' % ('\t'.join(print_names), attr_dict[key]))

    return outfile.close()


def build_species2prot2data_dict(species2prots, prot2data, species_token='|'):
    """ Build species-specific dictionaries that map species to prot name
    to data. """

    species2prot2data = {}
    for species in species2prots: 
        species2prot2data[species] = {}

    for prot in prot2data:
        species = prot.split(species_token)[0]
        species2prot2data[species][prot] = prot2data[prot]

    return species2prot2data
         

################################################################################
# AGING FUNCTIONS
################################################################################

def pull_out_unique_species(family_prots):
    """ Return a list of all species names that are present in the
    sequences names in the family_prots list. """

    species_list = [x.split('|')[0] for x in family_prots]

    return list(set(species_list))


def translate_tree(tree, leaf_dict):
    """ Replace the leaf labels in tree with whatever they map to in
    leaf_dict. """

    for leaf in tree.leaf_nodes():
        leaf_label = leaf.taxon.label.replace(' ', '_')
        leaf.taxon.label = leaf_dict[leaf_label]

    return tree


def run_count(species_tree, species2prots, family2prots, tree_fn,
              binary=False, species_token='|', program='AsymmetricWagner',
              params='-gain 1'):
    """ Return family2node2presense. """

    family2node2presence = {}

    species_order = [x.taxon.label for x in species_tree.leaf_nodes()]

    empty_species_counts = {}
    for species in species_order: empty_species_counts[species] = 0

    # write file with profile for each family for input to Count
    temp_fn = 'temp_%d.count_table' % random.randint(0, 100000)
    temp_file = open(temp_fn, 'w')
    temp_file.write('family\t%s\n' % '\t'.join(species_order))

    for family in family2prots:
        species_counts = empty_species_counts.copy()

        for prot in family2prots[family]:
            species = prot.split(species_token)[0]   
            if binary:
                species_counts[species] = 1
            else:
                species_counts[species] += 1

        profile = [species_counts[s] for s in species_order]

        temp_file.write('%s\t%s\n' %  (family, '\t'.join(['%d' % x for x in profile])))

    temp_file.close()

    # run Count
    program_name = 'ca.umontreal.iro.evolution.genecontent.%s %s' % (program, params)
    cmd = 'java -Xmx2048M -cp ~/opt/Count/Count.jar %s %s %s' % (program_name, 
                                                                 tree_fn, temp_fn)
    status, output = commands.getstatusoutput(cmd)

    # parse the output and build family2node2presence
    nodes = []
    for line in output.split('\n'):
        if "FAMILY" in line:
            t = line.split('\t')
            if t[1] == 'name':
                #nodes = [x.replace(' ', '_') for x in t[2:-4]]
                nodes = t[2:-4]
            else:
                family_name = t[1]
                counts = [int(x) for x in t[2:-4]]
                family2node2presence[family_name] = dict(zip(nodes, counts))
   
    os.remove(temp_fn)
        
    return family2node2presence


################################################################################
def age_proteins_count(species_tree, species2prots, family2prots, tree_fn,
                       species_token='|', recon_algo='wagner', gain_penalty=1,
                       exclude_species=[]):
    """ Return dictionaries mapping protein names to ages computed via
    the Count program. The program (and params) variable selects the
    algorithm used to reconstruct the family histories. There are two
    age dictionaries: one giving the depth of the LCA in the species
    tree and the other the label of the LCA node."""

    if recon_algo == 'wagner':
        family2node2presence = run_count(species_tree, species2prots, 
                                         family2prots, tree_fn, binary=True,
                                         program='AsymmetricWagner', 
                                         params='-gain %f' % gain_penalty)
    else:
        sys.exit('ERROR: %s is not a valid reconstruction algorithm.' % recon_algo)
        
    prot2age_label = {}
    prot2age_depth = {}

    # get depth of each species in tree
    max_depth = 0
    species_depth = {}
    for leaf in species_tree.leaf_nodes():
        leaf_depth = leaf.level()
        species_depth[leaf.taxon.label] = leaf_depth
        if leaf_depth > max_depth: max_depth = leaf_depth


    for family in family2prots:
        # compute age for each species in this family
        species2age = {}
        for leaf in species_tree.leaf_nodes():
            leaf_depth = leaf.level()
            leaf_name = str(leaf.taxon)

            if leaf_name in exclude_species: continue

            if leaf_name not in family2node2presence[family]:
                sys.exit('ERROR: %s not present Count results for %s.' % \
                             (leaf_name, family))

            family_profile = family2node2presence[family]
            if family_profile[leaf_name] == 0: continue

            # follow path to root until we hit a 0
            prev_node = (leaf_name, leaf_depth)
            for ancestor in leaf.ancestor_iter(inclusive=False):
                anc_name = ancestor.label
                anc_depth = ancestor.level()

                if anc_name not in family_profile:
                    sys.exit('ERROR: %s not present Count results for %s.' % \
                                 (anc_name, family))

                if family_profile[anc_name] == 1:
                    prev_node = (anc_name, anc_depth)
                else:
                    break

                #print '\t', ancestor.label, ancestor.taxon, ancestor.level()

            species2age[leaf_name] = prev_node


        for prot in family2prots[family]:
            species = prot.split(species_token)[0]
            if species in exclude_species: continue

            prot2age_depth[prot] = species_depth[species] - species2age[species][1]
            prot2age_label[prot] = species2age[species][0].replace(' ', '_')

    return prot2age_depth, prot2age_label

    
def age_proteins_dollo(species_tree, species2prots, family2prots, 
                        species_token='|'):
    """ Return dictionaries mapping protein names to ages computed via
    Dollo parsimony.  There are two age dictionaries: one giving the
    depth of the LCA in the species tree and the other the label of
    the LCA node."""

    prot2age_depth = {}
    prot2age_label = {}

    # get depth of each species in tree
    max_depth = 0
    species_depth = {}
    for leaf in species_tree.leaf_nodes():
        leaf_depth = leaf.level()
        species_depth[leaf.taxon.label] = leaf_depth
        if leaf_depth > max_depth: max_depth = leaf_depth

    # age each family
    for family in family2prots:
        # get all species that occur in the family
        taxon_labels = pull_out_unique_species(family2prots[family])

        # find the last common ancestor of these species
        print taxon_labels
        lca = species_tree.mrca(taxon_labels=taxon_labels)
        lca_label = lca.label
        if lca_label == None:
            if lca.taxon != None:
                lca_label = lca.taxon.label
            else:
                lca_label = '~'.join([x.taxon.label for x in lca.leaf_nodes()])

        lca_depth = lca.level()

        # assign ages to each family member
        # age is either the label of the lca or the distance from the
        # leaf to the lca
        for prot in family2prots[family]:
            species = prot.split(species_token)[0]

            prot2age_depth[prot] = species_depth[species] - lca_depth
            prot2age_label[prot] = lca_label.replace(' ', '_')

    return prot2age_depth, prot2age_label


def phylogenetic_profile(species_tree, species2prots, family2prots, 
                         binary=False, as_str=False, species_token='|'):
    """ Return dictionary mapping protein names to a list of counts in
    each species and a list giving the species ordering used in the
    profile.  If binary is True then return a 0/1 array, otherwise
    return counts of family members in other species. If as_str is
    True map prot to a string representation of the profile rather
    than a list of counts."""

    prot2profile = {}
    species_order = [x.taxon.label for x in species_tree.leaf_nodes()]

    empty_species_counts = {}
    for species in species_order:
        empty_species_counts[species] = 0

    # build profile for each family
    for family in family2prots:
        species_counts = empty_species_counts.copy()

        for prot in family2prots[family]:
            species = prot.split(species_token)[0]   
            if binary:
                species_counts[species] = 1
            else:
                species_counts[species] += 1

        profile = [species_counts[s] for s in species_order]
        
        for prot in family2prots[family]:
            if as_str:
                prot2profile[prot] = '\t'.join(['%d' % x for x in profile])
            else:
                prot2profile[prot] = profile

    return prot2profile, species_order


def parse_time_file(time_file):
    """ Return a dictionary mapping from labels to times (in MYA). """

    label2time = {}
    for line in open(time_file):
        if line[0] == '#': continue

        t = line[:-1].split()

        if len(t) != 2: continue

        label2time[t[0]] = float(t[1])

    return label2time


################################################################################

def main():

    # parse command line options and update variables
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ha:d:f:g:k:n:o:r:s:t:x:")
    except getopt.GetoptError:
        sys.exit(usage)

    if len(args) == 1:
        DATA_NAME = args[0].split('/')[-2]
        FAMILY_FILE = args[0] + DATA_NAME + '_families.txt'  
        TREE_FILE = args[0] + DATA_NAME + '.tree'
    elif len(args) == 2: 
        FAMILY_FILE = args[0]
        TREE_FILE = args[1]
    else:
        sys.exit(usage)

    RECONSTRUCTION_ALGO = 'dollo'  # -a
    RECONSTRUCTION_ALGOS = ['dollo', 'wagner']

    GAIN_PENALTY = 1 # -g

    OUTPUT_NAME = None # -n
    OUTPUT_DIR = None # -o

    OUTPUT_FORMATS = ['depth', 'label', 'time', 'profile'] # -f 
    POSSIBLE_OUTPUT_FORMATS = ['depth', 'label', 'time', 'profile']

    SPECIES_TOKEN = '|'  # -k, the character used to divide the species_id
                         # and protein_id in the protein name (e.g.,
                         # species_id|prot_id)

    SELECTED_SPECIES = None  # -s, only output ages for these species
                             # (comma separated list)

    TIME_FILE = None # -t, divergence times for nodes in the species
                     # tree in this file

    ADDITIONAL_ID_FILE = None  # -d, file mapping ids to aliases

    # -x, these are non-eukaryotic species in the DBs for which the trees
    # are not reliable and the analysis does not make sense
    EXCLUDE_SPECIES = ['SULSO', 'METAC', 'DEIRA', 'STRCO', 'THEMA', 'BRAJA', 
                       'ECOLI', 'PSEA7', 'GEOSL', 'AQUAE', 'CHLTA', 'GLOVI', 
                       'CHLAA', 'LEPIN', 'BACSU', 'BACTN']


    for opt, arg in opts:
        if opt == "-h":
            sys.exit('Use me properly.')

        elif opt == '-a':
            if arg in RECONSTRUCTION_ALGOS:
                RECONSTRUCTION_ALGO = arg
            else:
                sys.exit("ERROR: %s is not a valid reconstruction algorithm." % arg)

        elif opt == '-d':
            if os.path.isfile(arg):
                ADDITIONAL_ID_FILE = arg

        elif opt == '-f':
            in_formats = []
            for format in arg.split(','):
                if format in POSSIBLE_OUTPUT_FORMATS:
                    in_formats.append(format)
            if in_formats != []: OUTPUT_FORMATS = in_formats

        elif opt == '-g':
            try:
                GAIN_PENALTY = float(arg)
            except ValueError:
                sys.exit("ERROR: Gain penalty (%s) must be a positive real number." % arg)

        elif opt == '-k':
            if len(arg) != 1:
                sys.exit("ERROR: Token must be a single character. %s is too long." % arg)
            else:
                SPECIES_TOKEN = arg

        elif opt == '-n':
            OUTPUT_NAME = arg

        elif opt == '-o':
            OUTPUT_DIR = arg

        elif opt == '-s':
            SELECTED_SPECIES = arg.split(',')

        elif opt == '-t':
            if os.path.exists(arg):
                TIME_FILE = arg

        elif opt == '-x':
            EXCLUDE_SPECIES = arg.split(',')

    if OUTPUT_NAME == None:
        FAMILY_NAME = '.'.join(FAMILY_FILE.split('/')[-1].split('.')[:-1]).replace('_families', '')
        OUTPUT_NAME =  FAMILY_NAME + '_' + RECONSTRUCTION_ALGO 
        if RECONSTRUCTION_ALGO == 'wagner': OUTPUT_NAME += str(GAIN_PENALTY)

    if OUTPUT_DIR == None: OUTPUT_DIR = 'results/ages/%s/' % OUTPUT_NAME 

    if OUTPUT_DIR[-1] != '/': OUTPUT_DIR += '/'
    if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)

    PREAMBLE = '# %s %s\n#\n' % (' '.join(sys.argv), datetime.datetime.now())
    PREAMBLE += '# FAMILY_FILE = %s\n' % FAMILY_FILE
    PREAMBLE += '# TREE_FILE = %s\n' % TREE_FILE
    PREAMBLE += '# RECONSTRUCTION_ALGO = %s\n' % RECONSTRUCTION_ALGO
    PREAMBLE += '# GAIN_PENALTY = %.3g\n' % GAIN_PENALTY
    PREAMBLE += '# OUTPUT_NAME = %s\n' % OUTPUT_NAME
    PREAMBLE += '# OUTPUT_DIR = %s\n' % OUTPUT_DIR
    PREAMBLE += '# OUTPUT_FORMATS = %s\n' % OUTPUT_FORMATS
    PREAMBLE += '# TIME_FILE = %s\n' % TIME_FILE
    PREAMBLE += '# ADDITIONAL_ID_FILE = %s\n' % ADDITIONAL_ID_FILE
    PREAMBLE += '#\n'

    print PREAMBLE
    # Finished processing command line options

    # load data
    prot2family, family2prots, species2prots = \
        parse_families(FAMILY_FILE, species_token=SPECIES_TOKEN)

    species_tree = dendropy.Tree.get_from_path(TREE_FILE, schema="newick")

    # age proteins
    prot2age_depth = None; prot2age_label = None
    if RECONSTRUCTION_ALGO == 'dollo':
        prot2age_depth, prot2age_label = \
            age_proteins_dollo(species_tree, species2prots, family2prots,
                               species_token=SPECIES_TOKEN)

    elif RECONSTRUCTION_ALGO == 'wagner':
        prot2age_depth, prot2age_label = \
            age_proteins_count(species_tree, species2prots, family2prots, 
                               TREE_FILE, species_token=SPECIES_TOKEN,
                               recon_algo=RECONSTRUCTION_ALGO,
                               gain_penalty=GAIN_PENALTY,
                               exclude_species=EXCLUDE_SPECIES)

    else:
        sys.exit('ERROR: %s is not a valid reconstruction algorithm.' % \
                     RECONSTRUCTION_ALGO)

    # compute phylogenetic profiles
    prot2profile = None; profile_species_order = []
    if 'profile' in OUTPUT_FORMATS:
        prot2profile, profile_species_order = \
            phylogenetic_profile(species_tree, species2prots, family2prots, 
                                 as_str=True, species_token=SPECIES_TOKEN)

    # load additional aliases into dict
    id2alias = None
    if ADDITIONAL_ID_FILE:
        id2alias = {}
        for line in open(ADDITIONAL_ID_FILE):
            if line.startswith('#'): continue
            t = line[:-1].split('\t')
            #id2alias.setdefault(t[0], []).append(t[1])
            if t[0] in id2alias:
                id2alias[t[0]] += t[1:]
            else:
                id2alias[t[0]] = t[1:]


    # write files
    if 'label' in OUTPUT_FORMATS:
        # build species-specific dictionaries for convenience
        species2prot2age_label = \
            build_species2prot2data_dict(species2prots, prot2age_label, 
                                         species_token=SPECIES_TOKEN)

        for species in species2prot2age_label:
            if SELECTED_SPECIES and species not in SELECTED_SPECIES: continue
            if EXCLUDE_SPECIES and species in EXCLUDE_SPECIES: continue

            outfn = '%s%s_%s_age-label.protein_list' % (OUTPUT_DIR, species, OUTPUT_NAME)
            print 'Writing %s...' % outfn
            write_attribute_file_no_species(outfn, 
                                            species2prot2age_label[species], 
                                            preamble=PREAMBLE, 
                                            additional_id_dict=id2alias,
                                            species_token=SPECIES_TOKEN)

    if 'time' in OUTPUT_FORMATS and TIME_FILE:
        # build species-specific dictionaries for convenience
        species2prot2age_label = \
            build_species2prot2data_dict(species2prots, prot2age_label,
                                            species_token=SPECIES_TOKEN)
        label2time = parse_time_file(TIME_FILE)

        for species in species2prot2age_label:
            if SELECTED_SPECIES and species not in SELECTED_SPECIES: continue
            if EXCLUDE_SPECIES and species in EXCLUDE_SPECIES: continue

            label2time[species] = 0.0

            prot2age_time = {}
            for prot in species2prot2age_label[species]:
                label = species2prot2age_label[species][prot]
                if label in label2time:
                    prot2age_time[prot] = label2time[label]
                else:
                    s = 'WARNING: No time available in %s for label %s from %s.' % \
                        (TIME_FILE, label, prot)
                    print >> sys.stderr, s

            outfn = '%s%s_%s_age-time.protein_list' % (OUTPUT_DIR, species, OUTPUT_NAME)
            print 'Writing %s...' % outfn
            write_attribute_file_no_species(outfn, prot2age_time, 
                                            preamble=PREAMBLE,
                                            additional_id_dict=id2alias,
                                            species_token=SPECIES_TOKEN)


    if 'depth' in OUTPUT_FORMATS:
        # build species-specific dictionaries for convenience
        species2prot2age_depth = \
            build_species2prot2data_dict(species2prots, prot2age_depth,
                                            species_token=SPECIES_TOKEN)

        for species in species2prot2age_depth:
            if SELECTED_SPECIES and species not in SELECTED_SPECIES: continue
            if EXCLUDE_SPECIES and species in EXCLUDE_SPECIES: continue

            outfn = '%s%s_%s_age-depth.protein_list' % (OUTPUT_DIR, species, OUTPUT_NAME)
            print 'Writing %s...' % outfn
            write_attribute_file_no_species(outfn, species2prot2age_depth[species], 
                                            preamble=PREAMBLE, 
                                            additional_id_dict=id2alias,
                                            species_token=SPECIES_TOKEN)

    if 'profile' in OUTPUT_FORMATS:
        # build species-specific dictionaries for convenience
        species2prot2profile = \
            build_species2prot2data_dict(species2prots, prot2profile,
                                            species_token=SPECIES_TOKEN)

        for species in species2prot2profile:
            if SELECTED_SPECIES and species not in SELECTED_SPECIES: continue
            if EXCLUDE_SPECIES and species in EXCLUDE_SPECIES: continue

            outfn = '%s%s_%s_phylo-profile.protein_list' % (OUTPUT_DIR, species, 
                                                            OUTPUT_NAME)
            print 'Writing %s...' % outfn

            #add species labels to PREAMBLE comments
            profile_preamble = PREAMBLE
            profile_preamble += '#\n#protein\t%s\n' % ('\t'.join(profile_species_order))
            write_attribute_file_no_species(outfn, species2prot2profile[species], 
                                            preamble=profile_preamble, 
                                            additional_id_dict=id2alias,
                                            combine_ids='~',
                                            species_token=SPECIES_TOKEN)

    # End main()


if __name__ == "__main__": 
    main()

