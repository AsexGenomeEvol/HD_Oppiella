# -*- coding: utf-8 -*-

#######################################################################
### Patrick Tran Van : patrick.tranvan@gmail.com
###
### Removing contaminations after BlobTools. See the github for usage.
###
#######################################################################

import os
import re
import argparse
import sys
import glob
import os.path
import random
from Bio import SeqIO

#### HOW TO USE

## python contamination_filtration.py -s contamination_identification -i1 *.blobDB.table.txt
	
def contamination_identification(blob_table):
	"""
	Make a list of contaminant scaffolds.
	"""
	
	specie_contaminant_wide = []
	specie_contaminant_wide_dic = {}
            		
	blob_table_file = open(blob_table,"r")
	
	for entry in blob_table_file:
	
		if '#' in entry:	# skip header
			pass
		else:
			scaffold = re.split("[\r\t\n]", entry)[0]
			phylum = re.split("[\r\t\n]", entry)[9]
			phylum_hits = re.split("[\r\t\n]", entry)[12]
			# tax0=Streptophyta:837.0|Chordata:561.0|Eukaryota-undef:182.0|Ascomycota:165.0;

			genus = re.split("[\r\t\n]", entry)[21]
			specie = re.split("[\r\t\n]", entry)[25]
			
			not_contaminant_phylum_list = ["Annelida", "Arthropoda", "Brachiopoda", "Chaetognatha", "Chordata", "Cnidaria", "Echinodermata", "Hemichordata", "Mollusca", "Nematoda", "Nemertea", "no-hit", "Onychophora", "Placozoa", "Platyhelminthes", "Porifera", "Priapulida", "Rotifera", "Tardigrada", "Xenacoelomorpha"]	
							
			hit = False
			
			for not_contaminant_phylum in not_contaminant_phylum_list:
				if not_contaminant_phylum in phylum_hits:
					hit = True	# hits contain at least one metazoan or no-hit
					#print entry, phylum_hits
					break
					
			if hit == False:
				
				script = open("contaminant_scaffolds.txt","a")			
				script.write("{}\n".format(scaffold))
				script.close()
				
				script = open("contaminant_tax_list.txt","a")			
				script.write("{}\t{}\t{}\t{}\n".format(scaffold, phylum, genus, specie))
				script.close()
				
				specie_contaminant_wide.append(specie)
				specie_contaminant_wide_dic[specie] = phylum
			
	blob_table_file.close()

	for specie in list(set(specie_contaminant_wide)):
	
		script = open("unique_contaminant_specie_list.txt","a")			
		script.write("{}\n".format(specie))
		script.close()
		
		script = open("unique_contaminant_specie_phylum_list.txt","a")			
		script.write("{}\t{}\n".format(specie, specie_contaminant_wide_dic[specie]))
		script.close()
	
	print "\n\nTable has been written.\n"

													
def main(argv):
	
	mod=[]
	
	mod.append('\n%(prog)s -s contamination_identification -i1 <blob_table>')

	parser = argparse.ArgumentParser(prog = 'contamination_filtration.py',
                                 usage = "\n".join(mod))

	parser.add_argument('-s', action='store', dest='step_value',
	                    help='Step')
	                                                 	
	parser.add_argument('-i1', action='store', dest='input_value',
	                    help='Input 1')
	                    	                    	
	parser.add_argument('--version', action='version', version='%(prog)s 1.0')

	results = parser.parse_args()


	if results.step_value == "contamination_identification" and results.input_value:
		contamination_identification(results.input_value)
										
if __name__ == "__main__":
	main(sys.argv[1:])
