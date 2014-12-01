for thresh in *; do
	$(python /Users/virpatel/projects/vanderbilt-summer-2014/computing_age/ProteinHistorian/age_proteins.py -s hsa -t /Users/virpatel/projects/vanderbilt-summer-2014/computing_age/ProteinHistorian/example/time_tree_dates.txt /Users/virpatel/projects/vanderbilt-summer-2014/computing_age/phylo_tree_mirbase_mirviewer.txt $thresh)
done
