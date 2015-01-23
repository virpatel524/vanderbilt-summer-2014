echo $(python mirVIR.py  -age_file /Users/virpatel/projects/vanderbilt-summer-2014/main_script/mirbase_method2_ages.txt -gene_ages hgnc_names_2_age.txt -diseases microRNA_disease.txt -families /Users/virpatel/projects/vanderbilt-summer-2014/main_script/mirbase_parsed_families.txt -verified_targets /Users/virpatel/projects/vanderbilt-summer-2014/main_script/mirtarbase_verified_targets.txt)
# cd enrichment_results/

# for disease in /Users/virpatel/projects/vanderbilt-summer-2014/main_script/disease_mirnas/* ; do
# 	python /Users/virpatel/projects/vanderbilt-summer-2014/computing_age/ProteinHistorian/age_enrichment_analysis.py  -l /Users/virpatel/projects/vanderbilt-summer-2014/main_script/mirbase_method2_ages.txt /Users/virpatel/projects/vanderbilt-summer-2014/main_script/mirbase_method2_ages.txt $disease
# done




# python /Users/virpatel/projects/vanderbilt-summer-2014/computing_age/ProteinHistorian/age_enrichment_analysis.py  -l /Users/virpatel/projects/vanderbilt-summer-2014/main_script/mirbase_method2_ages.txt /Users/virpatel/projects/vanderbilt-summer-2014/main_script/mirbase_method2_ages.txt /Users/virpatel/projects/vanderbilt-summer-2014/main_script/mirnas_not_in_disease.txt
# python /Users/virpatel/projects/vanderbilt-summer-2014/computing_age/ProteinHistorian/age_enrichment_analysis.py  -l /Users/virpatel/projects/vanderbilt-summer-2014/main_script/mirbase_method2_ages.txt /Users/virpatel/projects/vanderbilt-summer-2014/main_script/mirbase_method2_ages.txt /Users/virpatel/projects/vanderbilt-summer-2014/main_script/mirnas_in_disease.txt
