for thresh in scrapes/*.txt; do
	echo $(python /Users/virpatel/projects/vanderbilt-summer-2014/computing_age/ProteinHistorian/age_proteins.py -s hsa -t /Users/virpatel/projects/vanderbilt-summer-2014/computing_age/ProteinHistorian/example/time_tree_dates.txt $thresh /Users/virpatel/textfiles/tree_manip/final_tree_replaced.txt )
done
