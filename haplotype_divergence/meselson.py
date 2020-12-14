# -*- coding: utf-8 -*-


#######################################################################
### Patrick Tran van : patrick.tranvan@gmail.com
###
### Extract alleles for Meselson study. See the github for usage.
### Version 02/11/20 - Subset version.
#######################################################################

import os
import re
import argparse
import sys
import glob
import os.path
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pprint
import collections
from itertools import combinations
import subprocess
import shutil
from pprint import pprint

#### HOW TO USE

## python $SC/meselson.py -s genome_af1_coverage_new -i1 subset.fasta -i2 subset.vcf -i3 subset.cov -i4 8 -o 2_As_b1v03_NAF.fasta
## python $SC/meselson.py -s fix_phaser -i1 As
## python $SC/meselson.py -s phase_v2 -i1 As_fix_input.txt -i2 100 -i3 3
## python $SC/meselson.py -s extract_allele_phaser_v2 -i1 P1I1_P1I2_P1I3_P2I1_P2I2_P2I3_P3I1_P3I2_P3I3.phase -i2 As_fix_input.txt -i3 2_As_b1v03_NAF.fasta 
## python $SC/meselson.py -s phase_v3 -i1 As_fix_input.txt -i2 100 -i3 combination.txt

			
def list_position_haplotype(variants, haplotype_A, haplotype_B):
	"""
	Output "Pos.A.B"
	"""
	
	list_pos_hap = []	
	
	variant = 0
	
	while variant < len(variants):
		list_pos_hap.append("{}.{}.{}".format(variants[variant].split(".")[1], haplotype_A[variant], haplotype_B[variant]))
		
		variant += 1
	
	return list_pos_hap
		

def genome_af1_coverage_new(genome_file, vcf_file, coverage_file, coverage_cutoff, output_file):
	"""
	Creating a new fasta file with new scaffolds that have the alternative bases for SNPs (1/1)
	and bases replaced by 'N' if the coverage < cut-off
	"""
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
	#coverage_cutoff = 10
		
	seq_record_list = []
	scaffold_dict = {}
            		
	vcfs = open(vcf_file,"r")
	
	# Parsing the VCF files
	
	for vcf in vcfs:

		if '#' in vcf:	# Skip the header
			pass
			
		else:
		
			scaffold = re.split("[\r\t\n]",vcf)[0]
			position = re.split("[\r\t\n]",vcf)[1] 
			ref_base = re.split("[\r\t\n]",vcf)[3]
			alt_base = re.split("[\r\t\n]",vcf)[4]
			af = re.split("[\r\t\n]",vcf)[9].split(":")[0]
			
			print "VCF: ", scaffold
			if len(ref_base) == 1 and len(alt_base) == 1 and af == "1/1":	# Fetch only SNP bi-allelic and alternative genotype 1/1 (the frequency that the alternate allele occurs in the genotypes)
			
				#print vcf
				
				if scaffold in scaffold_dict:
				
					scaffold_initial_seq = str(scaffold_dict[scaffold])
							
					# position - 1 because in python, coordinates start by 0.
							
					index = int(position)-1
						
					scaffold_new_seq = str(scaffold_initial_seq[:index]) + str(alt_base) + str(scaffold_initial_seq[index + 1:])
					scaffold_dict[scaffold] = scaffold_new_seq
				
				else:	# not existing yet, have to fetch from the genome.
						
					for seq_record in SeqIO.parse(genome_file, "fasta"):
						
						if seq_record.id == scaffold:	
							
							scaffold_initial_seq = str(seq_record.seq)	
																		
							index = int(position)-1
							
							scaffold_new_seq = str(scaffold_initial_seq[:index]) + str(alt_base) + str(scaffold_initial_seq[index + 1:])
							scaffold_dict[scaffold] = scaffold_new_seq
							
							
							break

					
	vcfs.close()
	
	# Add to the dictionary: scaffolds that have not been modified by VCF.
	
	for final_seq_record in SeqIO.parse(genome_file, "fasta"):
			
		if final_seq_record.id in scaffold_dict:
			pass
		else:
			scaffold_dict[final_seq_record.id] = str(final_seq_record.seq)
	
	#print scaffold_dict["1_Os_b1v03_scaf59714"]
	
	# Create the new fasta
	
	
	coverages = open(coverage_file,"r")
	
	for coverage in coverages:

		scaffold_coverage_file = re.split("[\r\t\n]",coverage)[0]
		print scaffold_coverage_file
		base_coverage_file = int(re.split("[\r\t\n]",coverage)[1]) - 1	# because coverage file is 1-based coordinates
		cov_coverage_file = int(re.split("[\r\t\n]",coverage)[2])
		
		if scaffold_coverage_file in scaffold_dict and cov_coverage_file < int(coverage_cutoff):

			scaffold_dict[scaffold_coverage_file] = str(scaffold_dict[scaffold_coverage_file][:base_coverage_file]) + 'N' + str(scaffold_dict[scaffold_coverage_file][base_coverage_file + 1:])
				
	coverages.close()
	
	# Order dictionary by key
	scaffold_dict_sort_by_key = collections.OrderedDict(sorted(scaffold_dict.items()))	
	
			
	for scf_order in scaffold_dict_sort_by_key:
		#print scf_order
		seq_record_list.append(SeqRecord(seq = Seq(scaffold_dict_sort_by_key[scf_order]), id=scf_order, name="", description=""))	
	#print scaffold_dict_sort_by_key
					
				
	SeqIO.write(seq_record_list, output_file, "fasta")
		
	print "\n\nFasta has been written.\n"


def fix_phaser(phaser_directory):

	"""
	Make new phaser files because some regions overlapped. 
	In that case: keep the longest region.
	"""
	
	fix_phaser_directory = os.path.basename(os.path.normpath(phaser_directory)) + "_fix"
	
	# Create new directory for new phaser files
	
	if os.path.exists(fix_phaser_directory):
		shutil.rmtree(fix_phaser_directory)
	
	os.makedirs(fix_phaser_directory)
		
	# Go through all phaser files
	
	data_repo = os.path.join(phaser_directory, "*haplotypic_counts*txt")
	
	data_file = 0
	while data_file < len(glob.glob(data_repo)):
	
		phaser_file = sorted(glob.glob(data_repo))[data_file]	
		sample_name = (phaser_file.split("/")[-1]).split(".")[0]

		output_file = os.path.join(fix_phaser_directory, sample_name + ".haplotypic_counts_mod.txt")

		regions = open(phaser_file,"r")
		
		phaser_dict = {}
		
		for region in regions:
	
			scaffold = re.split("[\r\t\n]",region)[0] 
			start = re.split("[\r\t\n]",region)[1] 
			stop = re.split("[\r\t\n]",region)[2]		
			rest = re.split("[\r\t\n]",region)[3:-1] 
			
			if scaffold == "contig":	# Skip the header
			
				script = open(output_file,"a")
				script.write("{}".format(region))
				script.close()
				
				pass
				
			else:
			
				phase_length = int(stop) - int(start) + 1	
				
				if scaffold in phaser_dict:
				
					overlap_longer = False
					finish = False
					overlap_exist = False
					
					while finish == False:
					
						scaffold_region = 0
						
						while scaffold_region < len(phaser_dict[scaffold]):	# Check if there are overlapped regions
										
							phase_present = range(int(phaser_dict[scaffold][scaffold_region][0]),int(phaser_dict[scaffold][scaffold_region][1])+1)
							phase_new = range(int(start),int(stop)+1)
							
							phase_present_set = set(phase_present)
							res_overlap = phase_present_set.intersection(phase_new)
							
							if len(res_overlap) > 0:	# Overlap exists
								
								overlap_exist = True
								
								if phaser_dict[scaffold][scaffold_region][2] > phase_length:	# Old region is longer
									
									overlap_longer = False
									scaffold_region += 1
									
								else:	# replace the region because old region is shorter
								
									overlap_longer = True
									#print phaser_dict[scaffold][scaffold_region]	# print the removed region			
									phaser_dict[scaffold].remove(phaser_dict[scaffold][scaffold_region])	# delete the shorter region
									
									break
							
							else:
								scaffold_region += 1
								
							#break
								
						
						if 	scaffold_region >= len(phaser_dict[scaffold]):
							finish = True
							
					if overlap_longer == True:
						phaser_dict[scaffold].append([start, stop, phase_length] + rest)	
						#print phaser_dict[scaffold]
						#print

					if overlap_exist == False:
						phaser_dict[scaffold].append([start, stop, phase_length] + rest)	
													
				else:
				
					phaser_dict[scaffold] = [[start, stop, phase_length] + rest]		
			
		regions.close()	
		
		# Writing the new phase files
		
		for key_scaffold, value_regions in phaser_dict.iteritems() :
			
			for value_region in value_regions:
			
				script = open(output_file,"a")
				script.write("{}\t{}\t{}\t{}\n".format(key_scaffold, value_region[0], value_region[1], ("\t").join(value_region[3:])))
				script.close()
				
		print sample_name + " ... done."
					
		data_file += 1
		

def bedops_command(comb, files_list):
	
	files_list_real = []
	
	for pop_indiv in comb:
		
		if pop_indiv == "P1I1":
			files_list_real.append(files_list[0])
			
		if pop_indiv == "P1I2":
			files_list_real.append(files_list[1])
			
		if pop_indiv == "P1I3":
			files_list_real.append(files_list[2])
			
		if pop_indiv == "P2I1":
			files_list_real.append(files_list[3])
			
		if pop_indiv == "P2I2":
			files_list_real.append(files_list[4])
			
		if pop_indiv == "P2I3":
			files_list_real.append(files_list[5])
			
		if pop_indiv == "P3I1":
			files_list_real.append(files_list[6])
			
		if pop_indiv == "P3I2":
			files_list_real.append(files_list[7])
			
		if pop_indiv == "P3I3":
			files_list_real.append(files_list[8])
	
	
	#cmd = '/scratch/beegfs/monthly/ptranvan/Software/bedops/2.4.35/bin/bedops --intersect {} > {}.intersect.bed'.format(" ".join(files_list_real), "_".join(comb))	
	cmd = 'bedops --intersect {}'.format(" ".join(files_list_real))	
	
	return cmd

def phase_output_v2(files_list, min_length, comb):

	files_list_real = []
	
	for pop_indiv in comb:
		
		if pop_indiv == "P1I1":
			files_list_real.append(files_list[0])
			
		if pop_indiv == "P1I2":
			files_list_real.append(files_list[1])
			
		if pop_indiv == "P1I3":
			files_list_real.append(files_list[2])
			
		if pop_indiv == "P2I1":
			files_list_real.append(files_list[3])
			
		if pop_indiv == "P2I2":
			files_list_real.append(files_list[4])
			
		if pop_indiv == "P2I3":
			files_list_real.append(files_list[5])
			
		if pop_indiv == "P3I1":
			files_list_real.append(files_list[6])
			
		if pop_indiv == "P3I2":
			files_list_real.append(files_list[7])
			
		if pop_indiv == "P3I3":
			files_list_real.append(files_list[8])
	
	#print files_list_real
	#print len(files_list_real)
	
	dict_phase = {}	# Create a dictionary of phased region to get the longest phase possible independently of individuals
	
	for phaser_file in files_list_real:
	
		regions = open(phaser_file,"r")
				
		for region in regions:
	
			scaffold = re.split("[\r\t\n]",region)[0] 
			
			if scaffold == "contig":	# Skip the header
				
				pass
			
			else:
			
				start = int(re.split("[\r\t\n]",region)[1])
				stop = int(re.split("[\r\t\n]",region)[2])
							
				if scaffold in dict_phase:
					dict_phase[scaffold]["indiv"].append(phaser_file)
					dict_phase[scaffold]["coord"].append(start)
					dict_phase[scaffold]["coord"].append(stop)
				
				else:
					dict_phase[scaffold]={}
					dict_phase[scaffold]["coord"]=[start, stop]
					dict_phase[scaffold]["indiv"]=[phaser_file]
					
		regions.close()
					
	#print dict_phase["2_As_b1v03_scaf02973"]
	
	# Order dictionary by key
	dict_phase_sort_by_key = collections.OrderedDict(sorted(dict_phase.items()))	

	# Writing phase file
	output_file = "_".join(comb) + ".phase"
	
	if os.path.isfile(output_file):
		os.remove(output_file)
	
		
	script = open(output_file,"a")
	script.write("Scaffold\tStart\tStop\tSize\n")
	script.close()
										
	for scf_order in dict_phase_sort_by_key:
	
		if len(set(dict_phase_sort_by_key[scf_order]["indiv"])) == len(files_list_real):	# all phaser files have this region

			max_phase = max(dict_phase_sort_by_key[scf_order]["coord"])
			min_phase = min(dict_phase_sort_by_key[scf_order]["coord"])
		
			if max_phase - min_phase >= int(min_length):

				script = open(output_file,"a")
				script.write("{}\t{}\t{}\t{}\n".format(scf_order, str(min_phase), str(max_phase), str(max_phase - min_phase)))
				script.close()
	
									
def subset_sample(phase_file_name, phaser_input):
			
	subset_sample_list = []
	files_list = []
	
	files = open(phaser_input,"r")
	
	# Parsing the list_files input
	
	for file in files:

		files_list.append(re.split("[\r\t\n]",file)[0])
							
	files.close()
	
	for sample in phase_file_name.split("_"):
		
		if sample == "P1I1":
			subset_sample_list.append(files_list[0])
			
		if sample == "P1I2":
			subset_sample_list.append(files_list[1])
			
		if sample == "P1I3":
			subset_sample_list.append(files_list[2])
			
		if sample == "P2I1":
			subset_sample_list.append(files_list[3])
			
		if sample == "P2I2":
			subset_sample_list.append(files_list[4])
			
		if sample == "P2I3":
			subset_sample_list.append(files_list[5])
			
		if sample == "P3I1":
			subset_sample_list.append(files_list[6])
			
		if sample == "P3I2":
			subset_sample_list.append(files_list[7])
			
		if sample == "P3I3":
			subset_sample_list.append(files_list[8])
			
	return subset_sample_list

		
def phase_v2(input_files, min_length, min_ind_pop):
	"""
	Output the longest regions where at least min_ind_pop overlap by min_length.
	"""
	
	files_list = []
	files = open(input_files,"r")
	
	# Parsing the list_files input
	
	for file in files:

		files_list.append(re.split("[\r\t\n]",file)[0])
							
	files.close()

	#print files_list
	
	files_list_var = ["P1I1", "P1I2", "P1I3", "P2I1", "P2I2", "P2I3", "P3I1", "P3I2", "P3I3"]
	
	# Get all combinations of files and length = min_ind_pop  
	combs = []
	
	if int(min_ind_pop) == 0:
	
		print "Processing ... intra population mode"
		
		comb = ["P1I1", "P1I2", "P1I3"]
		phase_output_v2(files_list, min_length, comb)
			
	elif int(min_ind_pop) == 1:
	
		ind = 3
		while ind < 10:
			combs += list(combinations(files_list_var, ind))
			ind += 1	 
	 
	elif int(min_ind_pop) == 2:
		
		ind = 6
		while ind < 10:
			combs += list(combinations(files_list_var, ind))
			ind += 1

	elif int(min_ind_pop) == 3:
		
		combs = combinations(files_list_var, 9) 
		
	else:
		print "BUG: min_ind_pop = 0, 1, 2 or 3"
		combs = []
			
	# Print the obtained combinations 
	for comb in list(combs): 

		if int(min_ind_pop) == 1:
		
			if len(filter(lambda x:'P1' in x, comb)) >= 1 and len(filter(lambda y:'P2' in y, comb)) >= 1 and len(filter(lambda z:'P3' in z, comb)) >= 1:
				
				print "Processing ... " + "_".join(comb)
				phase_output_v2(files_list, min_length, comb)

		if int(min_ind_pop) == 2:
		
			if len(filter(lambda x:'P1' in x, comb)) >= 2 and len(filter(lambda y:'P2' in y, comb)) >= 2 and len(filter(lambda z:'P3' in z, comb)) >= 2:
				
				print "Processing ... " + "_".join(comb)
				phase_output_v2(files_list, min_length, comb)
	    
		if int(min_ind_pop) == 3:
		
			if len(filter(lambda x:'P1' in x, comb)) == 3 and len(filter(lambda y:'P2' in y, comb)) == 3 and len(filter(lambda z:'P3' in z, comb)) == 3:
				print "Processing ... " + "_".join(comb)
				phase_output_v2(files_list, min_length, comb)
		
		#print comb

def iupac(base_A, base_B):
	"""
	Return IUPAC code for variantCount = 1
	https://www.bioinformatics.org/sms/iupac.html
	"""
	
	if (base_A == "A" and base_B == "G") or (base_A == "G" and base_B == "A"):
		return "R"

	if (base_A == "C" and base_B == "T") or (base_A == "T" and base_B == "C"):
		return "Y"
		
	if (base_A == "G" and base_B == "C") or (base_A == "C" and base_B == "G"):
		return "S"
		
	if (base_A == "A" and base_B == "T") or (base_A == "T" and base_B == "A"):
		return "W"
		
	if (base_A == "G" and base_B == "T") or (base_A == "T" and base_B == "G"):
		return "K"
		
	if (base_A == "A" and base_B == "C") or (base_A == "C" and base_B == "A"):
		return "M"
		
											
def extract_allele_phaser_v2(phase_file, phaser_input, genome_file):
	"""
	phase_file from phase_v2() function.
	Extract scaffolds and build haplotype A and B. for variantcount = 1: replace by iupac code.
	"""
	
	genome_prefix = (genome_file.split("/")[-1]).split(".")[0]
	
	phase_file_name = genome_prefix + "_" + (phase_file.split("/")[-1]).split(".")[0]
	
	if os.path.isfile("{}.fasta".format(phase_file_name)):
		os.remove("{}.fasta".format(phase_file_name))
		
	# Extracting the individuals 
	
	subset_sample_list = subset_sample(phase_file_name, phaser_input)
	
	#print subset_sample_list
	# ['As_fix/ASG1.haplotypic_counts_mod.txt', 'As_fix/ASG2.haplotypic_counts_mod.txt', 'As_fix/ASG7.haplotypic_counts.txt', 'As_fix/ASH1.haplotypic_counts_mod.txt', 'As_fix/ASH2.haplotypic_counts_mod.txt', 'As_fix/ASH3.haplotypic_counts_mod.txt', 'As_fix/ASS1.haplotypic_counts_mod.txt', 'As_fix/ASS2.haplotypic_counts_mod.txt', 'As_fix/ASS3.haplotypic_counts_mod.txt']

	# Create a dictionary of the genome
	scaffold_dict = {}
	
	for final_seq_record in SeqIO.parse(genome_file, "fasta"):

		scaffold_dict[final_seq_record.id] = str(final_seq_record.seq)
	
			
	seq_record_list = []
	
	
	# Go through phase file
	
	phase_regions = open(phase_file,"r")
	
	for phase_region in phase_regions:
	
		scaffold = re.split("[\r\t\n]",phase_region)[0]
		start = re.split("[\r\t\n]",phase_region)[1] 
		stop = re.split("[\r\t\n]",phase_region)[2]
		size = re.split("[\r\t\n]",phase_region)[3]
		
		
		haplotypic_counts_file = ""
		
		if scaffold == "Scaffold":	# Skip the header
			pass
		
		else:
			
			for haplotypic_counts_file in subset_sample_list:
				
				sample_name = (haplotypic_counts_file.split("/")[-1]).split(".")[0]
				#print sample_name
				
				haplotypic_counts = open(haplotypic_counts_file,"r")
				
				haplotype_A_seq = scaffold_dict[scaffold]				
				haplotype_B_seq = scaffold_dict[scaffold]
							
				# Parsing the haplotypic counts file (phaser output)
				
				for feature in haplotypic_counts:
			
					scaffold_haplotypic_counts = re.split("[\r\t\n]",feature)[0]
					start_haplotypic_counts = re.split("[\r\t\n]",feature)[1]
					stop_haplotypic_counts = re.split("[\r\t\n]",feature)[2]
					variants = re.split("[\r\t\n]",feature)[3] 
					variantCount = re.split("[\r\t\n]",feature)[4] 
					haplotype_A = re.split("[\r\t\n]",feature)[7]
					haplotype_B = re.split("[\r\t\n]",feature)[8]
					new_variants = []
					
						
					if scaffold_haplotypic_counts == "contig":	# Skip the header
						pass
						
					else:
						
						if scaffold == scaffold_haplotypic_counts:
											
							#print feature
							#print haplotype_A
							#print haplotype_B					
							
							new_variants = list_position_haplotype(variants.split(","), haplotype_A.split(","), haplotype_B.split(","))

							#print haplotypic_counts
							#print new_variants
							
							# Example: 2_As_b1v03_scaf00001.101435.A.G,2_As_b1v03_scaf00001.101466.G.A	2		0	G,G	A,A
							# become
							# Output "Pos.A.B"
							# ['101435.G.A', '101466.G.A']
							
							# for variantCount = 1
							# 2_As_b1v03_scaf00001	63276	63276	2_As_b1v03_scaf00001.63276.A.G	1		0	A	G
							# become 
							# ['63276.A.G']
						
							
							for variant_scaffold in new_variants:
					
								position = variant_scaffold.split(".")[0]
								haplotype_A_base = variant_scaffold.split(".")[1]	# haplotype A
								haplotype_B_base = variant_scaffold.split(".")[2]	# haplotype B
						
								#print position, haplotype_A_base, haplotype_B_base
						
								# position - 1 because in python, coordinates start by 0.
						
								index = int(position) - 1
								
								if int(variantCount) > 1:
						
									haplotype_A_seq = str(haplotype_A_seq[:index]) + str(haplotype_A_base) + str(haplotype_A_seq[index + 1:])
							
									haplotype_B_seq = str(haplotype_B_seq[:index]) + str(haplotype_B_base) + str(haplotype_B_seq[index + 1:])
								else:	# iupac code
							
									ambiguous_base = iupac(haplotype_A_base, haplotype_B_base)

									haplotype_A_seq = str(haplotype_A_seq[:index]) + ambiguous_base + str(haplotype_A_seq[index + 1:])
							
									haplotype_B_seq = str(haplotype_B_seq[:index]) + ambiguous_base + str(haplotype_B_seq[index + 1:])						
					
				haplotypic_counts.close()
				
				
				# Create the sequences
			
				# phase file is 1 based coordinates !
				
				haplotype_A_seq = str(haplotype_A_seq[int(start)-1:int(stop)])
				haplotype_B_seq = str(haplotype_B_seq[int(start)-1:int(stop)])
			
				id_seq = "{}_{}_{}_{}".format(scaffold, start, stop, sample_name)
			
				haplotype_A_allele = SeqRecord(seq = Seq(haplotype_A_seq), id=id_seq + "_haplotypeA", name="", description="")					
				haplotype_B_allele = SeqRecord(seq = Seq(haplotype_B_seq), id=id_seq + "_haplotypeB", name="", description="")	
			
				seq_record_list.append(haplotype_A_allele)
				seq_record_list.append(haplotype_B_allele)
			
				print id_seq
				
							
	phase_regions.close()	
	
	SeqIO.write(seq_record_list, "{}.fasta".format(phase_file_name), "fasta")	


def phase_v3(input_files, min_length, comb_file):
	"""
	Output the longest regions where at least x individuals overlap by min_length.
	input_files containt a list of individuals separated by '_' per line ex:
	
	P1I1_P1I2_P2I1
	P1I1_P1I2_P2I1_P2I2
	"""
	
	files_list = []
	files = open(input_files,"r")
	
	# Parsing the list_files input
	
	for file in files:

		files_list.append(re.split("[\r\t\n]",file)[0])
							
	files.close()

	# Parsing the comb_files input
	
	combinations = open(comb_file,"r")
	for combination in combinations:

		comb = re.split("[\r\t\n_]", combination)[:-1]
		print "Processing ... {}".format(comb)
		
		#print files_list, comb
		phase_output_v2(files_list, min_length, comb)
		

	combinations.close()
																					
def main(argv):
	
	mod=[]

	mod.append('\n%(prog)s -s genome_af1_coverage_new -i1 <genome_file> -i2 <vcf_file> -i3 <coverage_file> -i4 <coverage_cutoff> -o <output_file>')
	mod.append('%(prog)s -s fix_phaser -i1  <phaser_directory>')
	mod.append('%(prog)s -s genome_af1_coverage -i1 <genome_file> -i2 <vcf_file> -i3 <coverage_file> -o <output_file>')
	mod.append('%(prog)s -s phase_v2 -i1 <list_files> -i2 <min_length> -i3 <min_ind_pop>')
	mod.append('%(prog)s -s extract_allele_phaser_v2 -i1 <phase_file> -i2 <phaser_input> -i3 <genome>')
	mod.append('%(prog)s -s phase_v3 -i1 <list_files> -i2 <min_length> -i3 <comb_file>')
	
	parser = argparse.ArgumentParser(prog = 'meselson.py',
                                 usage = "\n".join(mod))

	parser.add_argument('-s', action='store', dest='step_value',
	                    help='Step')
	                                                 	
	parser.add_argument('-i1', action='store', dest='input_value',
	                    help='Input 1')

	parser.add_argument('-i2', action='store', dest='input2_value',
	                    help='Input 2')

	parser.add_argument('-i3', action='store', dest='input3_value',
	                    help='Input 3')

	parser.add_argument('-i4', action='store', dest='input4_value',
	                    help='Input 4')
	                    
	parser.add_argument('-e', action='store', dest='email_value',
	                    help='Email')	                
	                    	              	                    
	parser.add_argument('-o', action='store', dest='output_value',
	                    help='Output')
	                    	                    	                    	
	parser.add_argument('--version', action='version', version='%(prog)s 1.09')

	results = parser.parse_args()


	if results.step_value == "genome_af1_coverage_new" and results.input_value and results.input2_value and results.input3_value and results.input4_value and results.output_value:
		genome_af1_coverage_new(results.input_value, results.input2_value, results.input3_value, results.input4_value, results.output_value)

	if results.step_value == "fix_phaser" and results.input_value:
		fix_phaser(results.input_value)		
		
	if results.step_value == "phase_v2" and results.input_value and results.input2_value and results.input3_value:
		phase_v2(results.input_value, results.input2_value, results.input3_value)

	if results.step_value == "extract_allele_phaser_v2" and results.input_value and results.input2_value and results.input3_value:
		extract_allele_phaser_v2(results.input_value, results.input2_value, results.input3_value)

	if results.step_value == "phase_v3" and results.input_value and results.input2_value and results.input3_value:
		phase_v3(results.input_value, results.input2_value, results.input3_value)
	
																															
if __name__ == "__main__":
	main(sys.argv[1:])
