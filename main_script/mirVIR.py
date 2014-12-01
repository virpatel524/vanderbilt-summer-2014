import sys, os, string, numpy, matplotlib.pyplot as plt

mirna_ages = ''
gene_ages = ''
disease_associations = ''
family_associations = ''


args = sys.argv
for index, arg in enumerate(args):
	if arg == '-age_file':
		mirna_ages = args[index + 1]
		continue
	if arg == '-gene_ages':
		gene_ages = args[index + 1]
		continue
	if arg == '-diseases':
		continue
		disease_associations = args[index + 1]
	if 'tree' in arg:
		disease_associations = args[index + 1]


