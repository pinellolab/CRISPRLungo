#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pysam
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt
import math
import numpy as np
import pysam
import collections
from Bio import SeqIO
import sys

def analysis_function(control,edited,refernce, p_limit_value = 0.002,allowance_value = 0.05, pooling = True, Filter1 = False):
	# Path to SAM files
	samfile_path1 = edited
	samfile_path2 = control

	
	fasta_file = refernce

	def find_sub(read,fasta_file):
		reference_sequence = ""
		for record in SeqIO.parse(fasta_file, "fasta"):
			reference_sequence = str(record.seq)
			break

		# Open the SAM file
		substitution_list = []
	   
		read_seq = read.query_sequence
		
		# Iterate over the aligned portion of the read
		for read_idx, ref_idx in enumerate(read.get_reference_positions(full_length=True)):
			if ref_idx is None:
				continue
			read_base = read_seq[(read_idx)]
			if read_base != reference_sequence[ref_idx]:
				substitution_list.append(("Sub_"+str(reference_sequence[ref_idx])+">"+str(read_base),ref_idx,1))

		# Output the counts for each position
		return(substitution_list)

	def combine_deletions(tuples):
		# Function to check if two deletions should be combined
		if len(tuples)<=1:return(tuples)

		def should_combine(pos1, len1, pos2,len2):
			return (abs(pos2 - (pos1 + len1)) <= 0.05 * (len1+len2)/2 and len1>=5 and len2>=5)

		how_many_pop =0
		combined_tuples = []
		for i in range(0,len(tuples)):
			 tuple1 = tuples[i-how_many_pop]
			 if list(tuple1)[2]<=5: 
				  tuples.pop(i-how_many_pop)
				  how_many_pop+=1
				  combined_tuples.append(tuple1)
		i = 0
		while i < len(tuples)-1:
			if tuples[i][0] == "insertion" or tuples[i][0] == "deletion" and tuples[i+1][0] == tuples[i][0]:
				pos1, len1 = tuples[i][1], tuples[i][2]
				j = i + 1
				while j < (len(tuples)-1) and tuples[j][0] == tuples[i][0] and should_combine(pos1, len1, tuples[j][1], tuples[j][2]):
					len1 += tuples[j][2] + (tuples[j][1] - (pos1 + len1))
					j += 1
				combined_tuples.append((tuples[j][0], pos1, len1))
				i = j
			else:
				combined_tuples.append(tuples[i])
				i += 1
		return combined_tuples

	def quant_unique_indels(samfile_path):
		# Open the SAM file
		samfile = pysam.AlignmentFile(samfile_path, 'r')
		reference_sequence = ""
		for record in SeqIO.parse(fasta_file, "fasta"):
			reference_sequence = str(record.seq)
		# Initialize a dictionary to hold mutation information
		# Key: (type_of_mutation, position_of_mutation, length_of_mutation)
		# Value: count of occurrences
		mutations = {}
		#initilize dict where read number is key and value is list of mutations
		total_reads = 0
		counter = 0
		# Iterate over reads in the SAM file
		for read in samfile.fetch():
			mutations_in_read = []
			ref_pos = read.reference_start
			# Skip unmapped reads
			if read.is_unmapped:
				continue
			
			# Current position in the reference sequence
			
			if read.reference_start>100: continue
			if read.reference_end<(len(reference_sequence)-100): continue
			
			if read.mapping_quality <= 30:
				continue
	
			
			total_reads+=1
			# Iterate over CIGAR operations
			for operation, length in read.cigar:
				# Check for insertion (I) or deletion (D)
				if operation == 1:	# Insertion
					key = ('insertion', ref_pos, length)
					mutations_in_read.append(key)
				elif operation == 2 :  # Deletion
					key = ('deletion', ref_pos, length)
					mutations_in_read.append(key)
				# Update ref_pos based on operation
				if operation in [0, 2, 3]:	# Match/Mismatch, Deletion, N (Skipped region)
					ref_pos += length

				counter +=1
			#combined_mutations_in_read = combine_deletions(mutations_in_read)
			#if len(combined_mutations_in_read) != len(mutations_in_read):
			   # combined_mutations_in_read = combine_deletions(combined_mutations_in_read)
			
			#Substitution counter

			subs_in_read = find_sub(read,fasta_file)
			for sub in subs_in_read: mutations_in_read.append(sub)
			for mut in mutations_in_read:
				mutations[mut] = mutations.get(mut, 0) + 1

			
		# Close the SAM file
		samfile.close()
		return(mutations,total_reads)



	def custom_round(num,margin,length):
		if math.floor(num*length)  == 0: return(num)
		closest_value = round(num/math.floor(num*length))*math.floor(num*length)
		return int(closest_value)

	def pool(mut_dict, allowance = allowance_value):
		outputdict = {}
		for key in mut_dict.keys():
			if math.floor(list(key)[2]*allowance)  == 0: outputdict[key] = mut_dict[key]
			else: 
				newlength = custom_round(list(key)[2],allowance,list(key)[2])
				newpos = custom_round(list(key)[1],allowance,list(key)[2])
				type = list(key)[0]
				mutation = (type,newpos,newlength)
				outputdict[mutation] = outputdict.get(mutation, 0) + mut_dict[key]
		return(outputdict)

	def significant_mutations(control_dict, edited_dict, control_reads, edited_reads, p_limit = p_limit_value):
		"""
		Performs a Fisher's Exact Test on each mutation to compare its frequency between
		a control and an edited dataset. Returns the keys of the significant mutations.

		Parameters:
		- control_dict: Dictionary of mutations for the control dataset.
		- edited_dict: Dictionary of mutations for the edited dataset.
		- control_reads: Total number of reads in the control dataset.
		- edited_reads: Total number of reads in the edited dataset.

		Returns:
		- A list of keys for the mutations that are significantly different.
		"""
		

		significant_keys = []
		pvalues = []
		count_print = 0
		# Combine all keys from both dictionaries
		all_keys = set(control_dict.keys()) | set(edited_dict.keys())
		
		for key in all_keys:
			# Get the counts for this mutation in both datasets
			control_count = control_dict.get(key, 0)
			edited_count = edited_dict.get(key, 0)
			length = list(key)[2]
			if control_count/control_reads > edited_count/edited_reads:
				continue

			# Create the contingency table for this mutation
			table = [
				[control_count, control_reads - control_count],  # Control counts
				[edited_count, edited_reads - edited_count]		 # Edited counts
			]

			# Perform Fisher's Exact Test if control count is less than  or equal to 5 or chi2 if control_count >5
			if control_count <= 5:
				_, p_value = fisher_exact(table, alternative='two-sided')
			if control_count > 5:
				chi2, p_value, dof, expected = chi2_contingency(table)
			
		
			
			# Check for significance 
			if p_value<=p_limit:
					significant_keys.append(key)
					pvalues.append(p_value)
			
		return significant_keys,pvalues


	# Separate insertion and deletion counts for histogram plotting
	 

	def creat_dict_analysis(samfile_path, significant_keys):
		# Open the SAM file
		reference_sequence = ""
		for record in SeqIO.parse(fasta_file, "fasta"):
			reference_sequence = str(record.seq)
		samfile = pysam.AlignmentFile(samfile_path, 'r')

		# Initialize a dictionary to hold mutation information
		# Key: (type_of_mutation, position_of_mutation, length_of_mutation)
		# Value: count of occurrences
		mutations = {}
		#initilize dict where read number is key and value is list of mutations
		dict_of_reads = {}
		total_reads = 0
		set_sig_mut = set(significant_keys)
		# Iterate over reads in the SAM file
		reads_that_passed = []
		count_printed = 0
		for read in samfile.fetch():
			mutations_in_read = []
			# Skip unmapped reads
			if read.is_unmapped:
				continue
			#create more varible way to quantify end
			if read.reference_start>100 or read.reference_end<(len(reference_sequence)-100):
				continue

			if read.mapping_quality <= 30:
				continue
			reads_that_passed.append(str(read).split('	',1)[0])

			total_reads+=1
			#if total_reads%1000 == 0: print(read.reference_start)
			dict_of_reads[total_reads] = []
			# Current position in the reference sequence
			ref_pos = read.reference_start
			
			# Iterate over CIGAR operations
			for operation, length in read.cigar:
				# Check for insertion (I) or deletion (D)
				
				if operation == 1:	# Insertion
					key = ('insertion', ref_pos, length)

				elif operation == 2:  # Deletion
					key = ('deletion', ref_pos, length)

				if operation in [0, 2, 3]:	# Match/Mismatch, Deletion, N (Skipped region)
					ref_pos += length

				if operation not in [1,2]: continue

				mutations_in_read.append(key)
				
			subs_in_read = find_sub(read,fasta_file)
			for sub in subs_in_read: mutations_in_read.append(sub)
			
			#combined_mutations_in_read = combine_deletions(mutations_in_read)
			#if len(combined_mutations_in_read) != mutations_in_read: 
				#combined_mutations_in_read = combine_deletions(combined_mutations_in_read)
			if Filter1:
				if count_printed ==0:
					print("Proceeding With Statistical Tests")
					count_printed+=1
				for mut in mutations_in_read:
					
						allowance = allowance_value
						length1 = list(mut)[2]
						pos = list(mut)[1]
						type = list(mut)[0]
						if pooling:
							if (math.floor(length1*allowance)  == 0 or math.floor(list(mut)[1]*allowance)==0) :
								pos = list(mut)[1]
								type = list(mut)[0]
							else:
								length1 = custom_round(length1,allowance,list(mut)[2])
								pos = custom_round(list(mut)[1],allowance,list(mut)[2])
								type = list(mut)[0]
					
						mutation = (type,pos,length1)

						if mutation in set_sig_mut: dict_of_reads[total_reads].append(mut)
			else: 
				if count_printed ==0:
					print("Reporting All Mutations")
					count_printed+=1
				for mut in mutations_in_read:
					

					dict_of_reads[total_reads].append(mut)
				# Update ref_pos based on operation
			
				
		# Close the SAM file
		samfile.close()
		return(dict_of_reads,reads_that_passed)
	
	edited1,readedited1 = quant_unique_indels(samfile_path1)
	control1,readcontrol1 = quant_unique_indels(samfile_path2)

	if pooling: 
		control1 = pool(control1)
		edited1 = pool(edited1)

	significant_keys, pvalues = significant_mutations(control1, edited1, readcontrol1, readedited1)
	
	edited_dict_reads, reads_that_passed = creat_dict_analysis(samfile_path1,significant_keys)
	edited_dict_reads2, _ = creat_dict_analysis(samfile_path2,significant_keys)
	sub_count_edited = 0
	sub_count_control = 0
	for value in edited_dict_reads.values():
		for mutation in value:
			if mutation[0][0:3] == "Sub": sub_count_edited+=1
	for value in edited_dict_reads2.values():
		for mutation in value:
			if mutation[0][0:3] == "Sub": sub_count_control+=1
	print(sub_count_edited,sub_count_control)
	return(edited_dict_reads,edited_dict_reads2,reads_that_passed)
	
	


# In[3]:


import csv

def classify_mutations(mutations):
	if not mutations:
		return 'WT', 'None'
	
	def qualifier(type,length):
		if type == 'deletion':
			if length < 100:
				return 'Del'
			else:
				return 'Large_Del'
		elif type == 'insertion':
			if length < 20:
				return 'Ins'
			else:
				return 'Large_Ins'
		else: return(type)
	
	#2253:Sub_C>G {position}:Sub_{Ref_nt}>{Mut_nt}

	if len(mutations) > 1:
		return 'Complex', ','.join(f"{pos}_{pos+length-1}:{qualifier(type,length)}_{length}" for type, pos, length in mutations)
	
	type, pos, length = mutations[0]
	if type == 'deletion':
		if length < 100:
			return 'Del', f"{pos}_{pos+length-1}:Del_{length}"
		else:
			return 'Large_Del', f"{pos}_{pos+length-1}:Large_Del_{length}"
	elif type == 'insertion':
		if length < 20:
			return 'Ins', f"{pos}_{pos+length-1}:Ins_{length}"
		else:
			return 'Large_Ins', f"{pos}_{pos+length-1}:Large_Ins_{length}"
	else: return "Sub", f"{pos}_{pos+length-1}:{type}_{length}"
	

def process_mutations(mutations_dict, output_file,ids):
	with open(output_file, 'w', newline='') as file:
		writer = csv.writer(file, delimiter='\t')
		writer.writerow(['Read_id', 'Classification', 'Mutation_info'])
		
		for read_id, mutations in mutations_dict.items():
			read_id = ids[read_id-1]
			
			classification, mutation_info = classify_mutations(mutations)
			writer.writerow([read_id, classification, mutation_info])
