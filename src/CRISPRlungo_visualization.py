#!/usr/bin/env python
# coding: utf-8

import pysam, collections, sys, math, csv, re

import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import numpy as np

from Bio import SeqIO
from collections import Counter


def visualization_preprocess_regular(sam_file, fasta_file):

	def get_sequence_lengths(fasta_file):
		for record in SeqIO.parse(fasta_file, "fasta"):
			return len(record.seq)
	

	mutation_counts = collections.Counter()
	read_per_position = {}
	length = get_sequence_lengths(fasta_file)
	
	for pos in range(0, length):
		# Position in list corresponds to base A, T, G, C
		read_per_position[pos] = [0, 0, 0, 0]
	
	# Parse the reference sequence
	reference_sequence = ""
	for record in SeqIO.parse(fasta_file, "fasta"):
		reference_sequence = str(record.seq)
		break
	
	# Open the SAM file
	with pysam.AlignmentFile(sam_file, "r") as samfile:
		print_counter = 0
		for read in samfile:
			if not read.is_unmapped:
				if read.mapping_quality>=30:
					read_seq = read.query_sequence
					read_pos = read.reference_start
					if type(read_seq) != type('str'):
						print("continue")
						continue
					
					# Iterate over the aligned portion of the read
					adjuster = 0
					for read_idx, ref_idx in enumerate(read.get_reference_positions(full_length=True)):
						if ref_idx is None:
							continue
						read_base = read_seq[(read_idx)]

						# Increment the count for the respective base
						if read_base == 'A':
							read_per_position[ref_idx][0] += 1
						elif read_base == 'T':
							read_per_position[ref_idx][1] += 1
						elif read_base == 'G':
							read_per_position[ref_idx][2] += 1
						elif read_base == 'C':
							read_per_position[ref_idx][3] += 1
						else:
							print(ref_idx)
							break
	
	# Output the counts for each position

	return read_per_position

def align_count_plot(file_1, file_2, output_dir):

	f = open(file_1).readlines()
	labels = f[0].strip().split('\t')[:5]
	cnts = f[1].strip().split('\t')[:5]
	plotly_html = []

	x = ['Total_reads', 'Filter_unmapped', 'Filter_low_quality', 'Filter_short', 'Used_reads']
	y = cnts

	for i in range(5):
		y[i] = int(y[i])
		if 1 <= i <= 3:
			y[i] = y[i-1] - y[i]

	all_cnt = int(y[0])
	per = [round(int(i)/all_cnt*100,2) for i in y]

	fig = go.Figure()
	fig.add_trace(go.Bar(
		x=x,
		y=y,
		text=[f'{p:.2f}% ({int(yy):,})' for p, yy in zip(per, y)],
		textposition='outside',
		marker_color='grey',
		showlegend=False
	))
	
	# Customize the layout
	fig.update_layout(
		title='Alignment statistics',
		xaxis_title='Categories',
		yaxis_title='Number of Reads',
		yaxis=dict(tickformat=',', linecolor='black', linewidth=2),  # Black lines on the y-axis
		xaxis=dict(linecolor='black', linewidth=2),  # Black lines on the x-axis
		margin=dict(l=40, r=40, t=40, b=40),  # Adjust margins if needed
		plot_bgcolor='white'
	)
	
	# Show the plot
	plotly_html.append(fig.to_html(full_html=False, include_plotlyjs='cdn'))
	fig.write_image(f'{output_dir}/Control_alignment_statistics.png')
	fig.write_image(f'{output_dir}/Control_alignment_statistics.pdf')




	f = open(file_2).readlines()
	cnts = f[2].strip().split('\t')[:6]

	x = ['Total_reads', 'Filter_unmapped', 'Filter_low_quality', 'Filter_short', 'Filter_low_count', 'Used_reads']
	y = cnts

	for i in range(6):
		y[i] = int(y[i])
		if 1 <= i <= 4:
			y[i] = y[i-1] - y[i]

	all_cnt = int(y[0])
	per = [round(int(i)/all_cnt*100,2) for i in y]
	

	fig = go.Figure()
	fig.add_trace(go.Bar(
		x=x,
		y=y,
		text=[f'{p:.2f}% ({int(yy):,})' for p, yy in zip(per, y)],
		textposition='outside',
		marker_color='grey',
		showlegend=False
	))
	
	# Customize the layout
	fig.update_layout(
		title='Alignment statistics',
		xaxis_title='Categories',
		yaxis_title='Number of Reads',
		yaxis=dict(tickformat=',', linecolor='black', linewidth=2),  # Black lines on the y-axis
		xaxis=dict(linecolor='black', linewidth=2),  # Black lines on the x-axis
		margin=dict(l=40, r=40, t=40, b=40),  # Adjust margins if needed
		plot_bgcolor='white'
	)
	
	# Show the plot
	plotly_html.append(fig.to_html(full_html=False, include_plotlyjs='cdn'))
	fig.write_image(f'{output_dir}/Treated_alignment_statistics.png')
	fig.write_image(f'{output_dir}/Treated_alignment_statistics.pdf')


	return plotly_html

	

def regular_statistic_plot(tsv_file, fastq_file, output_dir):
	initial_reads = sum(1 for _ in SeqIO.parse(fastq_file, "fastq"))
	final_count = len(pd.read_csv(tsv_file, sep='\t', header=0))
	quantities = [initial_reads, initial_reads, final_count]  # Reads in input, Reads after preprocessing, Reads aligned
	labels = ["Reads in Input", "Reads after Preprocessing", "Reads Aligned"]
	
	# Calculate percentages
	total_reads = quantities[0]
	percentages = [(q / total_reads) * 100 for q in quantities]
	
	# Create the bar chart
	fig = go.Figure()
	
	# Add the bars
	fig.add_trace(go.Bar(
		x=labels,
		y=quantities,
		text=[f'{percentage:.2f}% ({quantity:,})' for percentage, quantity in zip(percentages, quantities)],  # Add percentage and count on top of the bars
		textposition='outside',
		marker_color='grey',
		showlegend=False
	))
	
	# Customize the layout
	fig.update_layout(
		title='Alignment statistics',
		xaxis_title='Categories',
		yaxis_title='Number of Reads',
		yaxis=dict(tickformat=',', linecolor='black', linewidth=2),  # Black lines on the y-axis
		xaxis=dict(linecolor='black', linewidth=2),  # Black lines on the x-axis
		margin=dict(l=40, r=40, t=40, b=40),  # Adjust margins if needed
		plot_bgcolor='white',
	)
	
	# Show the plot
	fig.write_image(f'{output_dir}/Alignment_statistics.png')
	fig.write_image(f'{output_dir}/Alignment_statistics.pdf')


def base_proportion(result, output_dir, reference, cv_pos, cv_pos_2, plot_window, show_all_between_allele):

	# Assume result is already defined and populated as in your initial script
	
	# Calculate total reads and base proportions
	
	window_st = cv_pos - plot_window
	window_ed = cv_pos + plot_window

	if cv_pos_2 != False:
		window_ed = cv_pos_2 + plot_window

	summed_results = {}
	base_proportions = {'A': [], 'T': [], 'G': [], 'C': []}
	
	for key in result.keys():
		if window_st > key or window_ed <= key:
			continue
		total_reads = sum(result[key])
		summed_results[key] = total_reads
		if total_reads > 0:
			base_proportions['A'].append(result[key][0] / total_reads)
			base_proportions['T'].append(result[key][1] / total_reads)
			base_proportions['G'].append(result[key][2] / total_reads)
			base_proportions['C'].append(result[key][3] / total_reads)
		else:
			base_proportions['A'].append(0)
			base_proportions['T'].append(0)
			base_proportions['G'].append(0)
			base_proportions['C'].append(0)
	
	# Base colors
	base_colors = {'A': 'honeydew',
		'T': 'mistyrose',
		'G': 'lightyellow',
		'C': 'aliceblue'}
	
	ref_seq_dictionary = {'A': [], 'T': [], 'G': [], 'C': []}
	
	poscounter =0 
	for nt in str(reference):
		if nt not in 'ATGC':
			continue
		ref_seq_dictionary[nt].append(poscounter)
		poscounter +=1
	
	# Create figure
	fig = go.Figure()
	
	# Add total reads bar
	fig.add_trace(go.Bar(
		x=list(summed_results.keys()), 
		y=list(summed_results.values()), 
		marker_color='black', 
		marker_line_width=0, 
		name='Total Reads'
	))
	
	ajustment = []
	for i in range(0,len(base_proportions["A"])):
		ajustment.append(0)
	maxreads = max(summed_results.values())

	# Add base proportion bars
	for base in ['A', 'T', 'G', 'C']:
		fig.add_trace(go.Bar(
			x=np.add(list(summed_results.keys()),-0.5), 
			y=[p * maxreads for	p in base_proportions[base]], 
			base=[(-maxreads)+adj for adj in ajustment],
			marker_color=base_colors[base], 
			marker_line_width=0, 
			name=f'{base} Proportion',
			width=0.8, 
			offset=0.1,  
		))
		ajustment=[adj+p * maxreads for adj,p in zip(ajustment,base_proportions[base])]
	
	
	"""for base in ['A', 'T', 'G', 'C']:
		fig.add_trace(go.Bar(
			x=np.add(ref_seq_dictionary[base], -0.5),
			y=[maxreads/10 for a in ref_seq_dictionary[base]], 
			base=[(-maxreads/4) for a in ref_seq_dictionary[base]],
			marker_color=base_colors[base], 
			showlegend=False,
			marker_line_width=0, 
			width=0.8, 
			offset=0.1,  
		))"""
		
	fig.update_layout(
		title='Reads by Position with Base Proportions',
		xaxis_title='Position',
		yaxis_title='Reads',
		barmode='stack',  # Stack bars to show proportions
		dragmode='pan',  # Allows dragging to pan around
		xaxis_rangeslider_visible=True,	# Adds a range slider at the bottom
		height = 400,
		plot_bgcolor='white',    
    	paper_bgcolor='white'    
	)
	
	# Show the figure
	fig.write_image(f'{output_dir}/Base_proportions.png')
	fig.write_image(f'{output_dir}/Base_proportions.pdf')
	return fig.to_html(full_html=False, include_plotlyjs='cdn')


def regular_accuracy_plot(reference, result, output_dir):
	
	x_vals = []
	correct_score = []
	for ref_idx in range(0,len(reference)):
		if sum(result[ref_idx]) == 0:
			continue
		x_vals.append(ref_idx)
		read_base = reference[ref_idx]
		if read_base == 'A':
			correct_score.append(result[ref_idx][0]/sum(result[ref_idx]))
		elif read_base == 'T':
			correct_score.append(result[ref_idx][1]/sum(result[ref_idx]))
		elif read_base == 'G':
			correct_score.append(result[ref_idx][2]/sum(result[ref_idx]))
		elif read_base == 'C':
			correct_score.append(result[ref_idx][3]/sum(result[ref_idx]))
	
	plt.scatter(x_vals, correct_score, color='red')
	theta = math.pi
	for i in range(len(x_vals)):
		if correct_score[i] < 0.5:
			plt.annotate(
				str(x_vals[i]), 
				(x_vals[i], correct_score[i]), 
				textcoords="offset points", 
				xytext=(45*math.cos(theta), 45*math.sin(theta)),  # Offset position
				ha='center',
				arrowprops=dict(
					arrowstyle="-", 
					color='black',
					lw=0.5	# Line width
				)
			)
			theta-=math.pi/6
	plt.title('Position vs accuracy')
	plt.xlabel('Position')
	plt.ylabel('Proportion of bases at position that match refeernce sequence')
	plt.grid(True)
	plt.savefig(f'{output_dir}/Regular_accuracy.png')
	plt.savefig(f'{output_dir}/Regular_accuracy.pdf')




def mutation_pie_chart(tsv_file, output_dir):
	def update_allele_dict(mutations, percent_of_alleles_dict):
		# Keep track of whether Del/LargeDel and Ins/LargeIns have been counted
		del_counted = False
		ins_counted = False
		mutations = mutations.split(",")
		if mutations[0]=="None": return(percent_of_alleles_dict)
		#print(mutation)
		for mutation in mutations:
			# Split the mutation into its components
			#print(mutations)
			mutation_type = mutation.split(':')[1]	# Get the mutation type, e.g., 'Sub_A->G_1'
			
			if 'Del' in mutation_type:
				if 'LargeDel' in mutation_type and not del_counted:
					percent_of_alleles_dict['LargeDel'] += 1
					del_counted = True
				elif 'Del' in mutation_type and not del_counted:
					percent_of_alleles_dict['Del'] += 1
					del_counted = True	# Ensure Del and LargeDel are not double counted
			elif 'Ins' in mutation_type:
				if 'LargeIns' in mutation_type and not ins_counted:
					percent_of_alleles_dict['LargeIns'] += 1
					ins_counted = True
				elif 'Ins' in mutation_type and not ins_counted:
					percent_of_alleles_dict['Ins'] += 1
					ins_counted = True	# Ensure Ins and LargeIns are not double counted
			elif 'Sub' in mutation_type:
				percent_of_alleles_dict['Sub'] += 1
			elif 'WT' in mutation_type:
				percent_of_alleles_dict['WT'] += 1

		return percent_of_alleles_dict

	# Initialize a list to store values from the second column
	second_column_values = []
	
	# Read the TSV file
	percent_of_alleles_dict = {"WT":0,'Del':0,'Ins':0,'LargeDel':0,"LargeIns":0,"Sub":0,'Inv':0, 'Complex': 0, 'ComplexWithLargeMut': 0}
	pie_colors = ['rgb(179,179,179)', '#DC3912', '#636EFA', '#FF9900', '#19D3F3', 'rgb(102,102,102)', '#00CC96', '#222A2A', '#6D6D6D']
	with open(tsv_file, mode='r') as file:
		next(file)  # Skip the header row
		for row in file:
			row = row.strip().split('\t')
			value = row[0]
			percent_of_alleles_dict[row[0]] += int(row[1])
			#percent_of_alleles_dict = update_allele_dict(row[2],percent_of_alleles_dict)
			try:
				# Attempt to convert to float and round if possible
				value = round(float(value), 1)
			except ValueError:
				# If conversion fails, keep the value as is (likely a string like 'WT')
				pass
			second_column_values.append(value)

	# Count occurrences of each value
	value_counts = Counter(second_column_values)
	
	# Dictionary for the first pie chart (WT vs Others)
	wt_count = percent_of_alleles_dict['WT']
	other_count = sum(percent_of_alleles_dict.values()) - wt_count
	wt_vs_other = {'UNMODIFIED': wt_count, 'MODIFIED': other_count}
	
	# Dictionary for the second pie chart (each key in the dictionary)
	detailed_counts = dict(value_counts)
	
	# Create the first pie chart (WT vs Other)
	fig1 = go.Figure(data=[go.Pie(labels=list(wt_vs_other.keys()), values=list(wt_vs_other.values()))])
	fig1.update_layout(showlegend=False, title='Proportion of WT vs Other')
	
	fig1.update_traces(
			textposition="outside",
			texttemplate="%{label} %{value} %{percent:.2%}",
	)
	
	percent_of_alleles_dict["WT"] = wt_count
	# Create the second pie chart (detailed counts)
	fig2 = go.Figure(data=[go.Pie(labels=list(percent_of_alleles_dict.keys()), 
		values=list(percent_of_alleles_dict.values()),  marker=dict(colors=pie_colors), domain={'x': [0, 0.7], 'y': [0, 0.7]})])
	fig2.update_layout(
		showlegend=False,
		title='Detailed Counts',
		legend=dict(
			x=1.1,
			y=0.5
		)
	)
	
	fig2.update_traces(
			textposition="outside",
			texttemplate="%{label} %{value} %{percent:.2%}",
			textfont=dict(size=16)
	)

	fig3 = go.Figure(data=[go.Pie(labels=list(percent_of_alleles_dict.keys()), values=list(percent_of_alleles_dict.values()))])
	fig3.update_layout(showlegend=False, title='Frequency of Alleles')

	fig3.update_traces(
			textposition="outside",
			texttemplate="%{label} %{value} %{percent:.2%}",
	)

	# Display the plots
	html_mut_pie = [fig1.to_html(full_html=False, include_plotlyjs='cdn')]
	fig1.write_image(f'{output_dir}/Mutation_pie_chart.png')
	fig1.write_image(f'{output_dir}/Mutation_pie_chart.pdf')
	html_pat_pie = [fig2.to_html(full_html=False, include_plotlyjs='cdn')]
	fig2.write_image(f'{output_dir}/Pattern_pie_chart.png')
	fig2.write_image(f'{output_dir}/Pattern_pie_chart.pdf')
	html_allele_pie = [fig3.to_html(full_html=False, include_plotlyjs='cdn')]
	fig3.write_image(f'{output_dir}/Percent_of_alleles_pie_chart.png')
	fig3.write_image(f'{output_dir}/Percent_of_alleles_pie_chart.pdf')

	return html_mut_pie, html_pat_pie, html_allele_pie



def custom_mutation_pie_chart(cnt_dict, output_dir):

	value_list = []
	key_list = []
	#for i in sorted(cnt_dict.items(), key = lambda x: x[1], reverse=True):
	for i in cnt_dict.items():
		value_list.append(i[1])
		key_list.append(i[0])

	fig1 = go.Figure(data=[go.Pie(labels=key_list, values=value_list, sort=False)])
	fig1.update_layout(showlegend=False, title='Proportion of Customized Mutations')
	
	fig1.update_traces(
			textposition="outside",
			texttemplate="%{label} %{value} %{percent:.2%}",
			
	)

	# Display the plots
	html_mut_pie = [fig1.to_html(full_html=False, include_plotlyjs='cdn')]
	fig1.write_image(f'{output_dir}/Custom_mutation_pie_chart.png')
	fig1.write_image(f'{output_dir}/Custom_mutation_pie_chart.pdf')

	return html_mut_pie


# Small/Large Insertion pos graph

def indel_per_position(tsv_file, reference, output_dir):

	
	df = pd.read_csv(tsv_file, sep='\t', header=None, names=['Classification', 'Count', 'Percentage', 'Insert_info', 'Induced_info', 'Mutation_info'])
	
	# Function to parse mutation info
	deletion_PerPos_dictionary = {}
	insertions_PerPos_dictionary = {}
	substitutions_PerPos_dictionary = {}
	
	# Initialize dictionaries
	for i in range(0, len(reference)): 
		deletion_PerPos_dictionary[i] = 0
		insertions_PerPos_dictionary[i] = 0  
		substitutions_PerPos_dictionary[i] = 0
	total_reads = 0
	# Extract and count mutations
	for idx, row in df.iterrows():
		if row['Count'] == 'Count':
			continue
		
		total_reads += int(row['Count'])
		
		if row['Classification'] != 'WT':
			parsed = row['Mutation_info'].split(",")
			for mut in parsed:
				if "Del" in str(mut):
					start_pos = int(re.match(r"(\d+)_", mut).group(1))
					end_pos = int(re.search(r"_(\d+)", mut).group(1))
					for position in range(start_pos, end_pos+1):
						deletion_PerPos_dictionary[position] += int(row['Count'])
				if "Ins" in str(mut):
					start_pos = int(re.match(r"(\d+)_", mut).group(1))
					end_pos = int(re.search(r"_(\d+)", mut).group(1))
					for position in range(start_pos, end_pos+1):
						insertions_PerPos_dictionary[position] += int(row['Count'])
				if "Sub" in str(mut):
					start_pos = int(re.match(r"(\d+)_", mut).group(1))
					end_pos = int(re.search(r"_(\d+)", mut).group(1))
					for position in range(start_pos, end_pos+1):
						substitutions_PerPos_dictionary[position] += int(row['Count'])
	
	max_deletions = max(deletion_PerPos_dictionary.values())
	max_insertions = max(insertions_PerPos_dictionary.values())
	max_sub = max(substitutions_PerPos_dictionary.values())
	max_value = max(max_deletions, max_insertions,max_sub)
	
	# Define Y-axis ticks
	y_ticks = [i * max_value / 5 for i in range(6)]
	y_tick_labels = [f"{(tick / total_reads) * 100:.2f}% ({int(tick)})" for tick in y_ticks]
	
	# Create the plot
	fig = go.Figure()
	
	# Add the smooth line for deletions
	fig.add_trace(go.Scatter(x=list(deletion_PerPos_dictionary.keys()), 
							 y=list(deletion_PerPos_dictionary.values()), 
							 mode='lines', 
							 line=dict(color='#8172b3'),
							 name='Deletions'))
	
	# Add the smooth line for insertions
	fig.add_trace(go.Scatter(x=list(insertions_PerPos_dictionary.keys()), 
							 y=list(insertions_PerPos_dictionary.values()), 
							 mode='lines', 
							 line=dict(color='#c44f51'),  # Different color for insertions
							 name='Insertions'))
	
	# Add Subs
	fig.add_trace(go.Scatter(x=list(substitutions_PerPos_dictionary.keys()), 
							 y=list(substitutions_PerPos_dictionary.values()), 
							 mode='lines', 
							 line=dict(color='green'),	# Different color for insertions
							 name='Substitutions'))
	
	# Customize the layout
	fig.update_layout(
		title='Deletions and Insertions per Position',
		xaxis_title='Position',
		yaxis_title='Sequence %',
		showlegend=True,
		dragmode='pan',  # Allows dragging to pan around
		xaxis_rangeslider_visible=True	# Adds a range slider at the bottom
	)
	
	fig.update_yaxes(tickvals=y_ticks, ticktext=y_tick_labels)
	
	# Show the plot
	fig.write_image(f'{output_dir}/Deletion_and_insertions_per_position.png')
	fig.write_image(f'{output_dir}/Deletion_and_insertions_per_position.pdf')
	return fig.to_html(full_html=False, include_plotlyjs='cdn')


# Insertion average length vs position

# In[32]:


def Insertion_length(tsv_file, output_dir):
	
	# Read the data
	df = pd.read_csv(tsv_file, sep='\t', header=None, names=['Classification', 'Count', 'Percentage', 'Insert_info', 'Induced_info', 'Mutation_info'])
	
	# Initialize dictionaries
	insertions_PerPos_dictionary = {}
	insertion_lengths = {}
	
	# Extract and count mutations
	for idx, row in df.iterrows():
		if row['Classification'] != 'WT':
			parsed = row['Mutation_info'].split(",")
			for mut in parsed:
				if "Ins" in str(mut):
					start_pos = int(re.match(r"(\d+)_", mut).group(1))
					length = int(re.search(r"Ins_(\d+)", mut).group(1))
					if start_pos not in insertions_PerPos_dictionary:
						insertions_PerPos_dictionary[start_pos] = 0
						insertion_lengths[start_pos] = []
					insertions_PerPos_dictionary[start_pos] += int(row['Count'])
					insertion_lengths[start_pos].append(length)
	
	# Calculate average insertion length at each position
	average_insertion_lengths = {pos: (sum(lengths) / len(lengths)) if lengths else 0 for pos, lengths in insertion_lengths.items()}
	
	# Create the plot
	fig = go.Figure()
	
	# Add the histogram bars for average insertion length
	fig.add_trace(go.Bar(
		x=list(average_insertion_lengths.keys()), 
		y=list(average_insertion_lengths.values()), 
		marker=dict(color='red'),
		width=0.5,
		opacity=0.7,
		name='Average Insertion Length',
		hoverinfo='x+y',
		showlegend=False
	))
	
	# Add the small red circles on top of the bars
	fig.add_trace(go.Scatter(
		x=list(average_insertion_lengths.keys()), 
		y=list(average_insertion_lengths.values()), 
		mode='markers', 
		marker=dict(color='red', size=5),
		name='Average Insertion Length Markers',
		showlegend=False
	))
	
	# Customize the layout
	fig.update_layout(
		title='Average Insertion Length at Each Position',
		xaxis_title='Position',
		yaxis_title='Average Insertion Length',
		showlegend=True,
		bargap=0.1,  # Adjust the gap between bars
		margin=dict(l=0, r=0, t=50, b=50),	# Adjust margins if needed
		xaxis_rangeslider_visible=True	# Adds a range slider at the bottom
	)
	
	# Show the plot
	fig.write_image(f'{output_dir}/insertion_length.png')
	fig.write_image(f'{output_dir}/insertion_length.pdf')
	return fig.to_html(full_html=False, include_plotlyjs='cdn')


def Deletion_length(tsv_file, output_dir):
	
	# Read the data
	df = pd.read_csv(tsv_file, sep='\t', header=None, names=['Classification', 'Count', 'Percentage', 'Insert_info', 'Induced_info', 'Mutation_info'])
	
	# Initialize dictionaries
	deletions_PerPos_dictionary = {}
	deletion_lengths = {}
	
	# Extract and count mutations
	for idx, row in df.iterrows():
		if row['Classification'] != 'WT':
			parsed = row['Mutation_info'].split(",")
			for mut in parsed:
				if "Del" in str(mut):
					start_pos = int(re.match(r"(\d+)_", mut).group(1))
					end_pos = int(re.search(r"_(\d+)", mut).group(1))
					length = end_pos - start_pos + 1
					if start_pos not in deletions_PerPos_dictionary:
						deletions_PerPos_dictionary[start_pos] = 0
						deletion_lengths[start_pos] = []
					deletions_PerPos_dictionary[start_pos] += int(row['Count'])
					deletion_lengths[start_pos].append(length)
	
	# Calculate average deletion length at each position
	average_deletion_lengths = {pos: (sum(lengths) / len(lengths)) if lengths else 0 for pos, lengths in deletion_lengths.items()}
	
	# Create the plot
	fig = go.Figure()
	
	# Add the histogram bars for average deletion length
	fig.add_trace(go.Bar(
		x=list(average_deletion_lengths.keys()), 
		y=list(average_deletion_lengths.values()), 
		marker=dict(color='blue'),
		width=0.5,
		opacity=0.7,
		name='Average Deletion Length',
		hoverinfo='x+y',
		showlegend=False
	))
	
	# Add the small blue circles on top of the bars
	fig.add_trace(go.Scatter(
		x=list(average_deletion_lengths.keys()), 
		y=list(average_deletion_lengths.values()), 
		mode='markers', 
		marker=dict(color='blue', size=5),
		name='Average Deletion Length Markers',
		showlegend=False
	))
	
	# Customize the layout
	fig.update_layout(
		title='Average Deletion Length at Each Position',
		xaxis_title='Position',
		yaxis_title='Average Deletion Length',
		showlegend=True,
		bargap=0.1,  # Adjust the gap between bars
		margin=dict(l=0, r=0, t=50, b=50),	# Adjust margins if needed
		xaxis_rangeslider_visible=True	# Adds a range slider at the bottom
	
	)
	
	# Show the plot
	fig.write_image(f'{output_dir}/deletion_length.png')
	fig.write_image(f'{output_dir}/deletion_length.pdf')
	return fig.to_html(full_html=False, include_plotlyjs='cdn')

def Deletion_count_length(tsv_file, output_dir):

	# Read the data
	df = pd.read_csv(tsv_file, sep='\t', header=None, names=['Classification', 'Count', 'Percentage', 'Insert_info', 'Induced_info', 'Mutation_info'])
	
	# Initialize dictionaries
	deletion_lengths_count = {}
	
	# Extract and count mutations
	for idx, row in df.iterrows():
		if row['Classification'] != 'WT':
			parsed = row['Mutation_info'].split(",")
			for mut in parsed:
				if "Del" in str(mut):
					start_pos = int(re.match(r"(\d+)_", mut).group(1))
					end_pos = int(re.search(r"_(\d+)", mut).group(1))
					length = end_pos - start_pos + int(row['Count'])
					if length not in deletion_lengths_count:
						deletion_lengths_count[length] = 0
					deletion_lengths_count[length] += int(row['Count'])
	
	# Create the plot
	fig = go.Figure()
	
	# Add the bars for deletion length counts
	fig.add_trace(go.Bar(
		x=list(deletion_lengths_count.keys()), 
		y=list(deletion_lengths_count.values()), 
		marker=dict(color='blue'),
		width=0.8,
		opacity=1,
		name='Deletion Length Count',
		hoverinfo='x+y',
		showlegend=False
	))
	
	# Customize the layout
	fig.update_layout(
		title='Count of Deletions vs. Deletion Length',
		xaxis_title='Deletion Length',
		yaxis_title='Count of Deletions (Log Scale)',
		showlegend=True,
		yaxis_type='log',  # Set y-axis to logarithmic scale
		bargap=0.1,  # Adjust the gap between bars
		margin=dict(l=0, r=0, t=50, b=50),	# Adjust margins if needed
		xaxis_rangeslider_visible=True
	)
	
	# Show the plot
	fig.write_image(f'{output_dir}/Deletion_count_length.png')
	fig.write_image(f'{output_dir}/Deletion_count_length.pdf')
	return fig.to_html(full_html=False, include_plotlyjs='cdn')


def insertion_count_length(tsv_file, output_dir):

	# Example data
	tsv_file = tsv_file
	
	# Read the data
	df = pd.read_csv(tsv_file, sep='\t', header=None, names=['Read_id', 'Classification', 'Mutation_info', 'Insert_info', 'Induced_info'])
	
	# Initialize dictionaries
	insertion_lengths_count = {}
	
	# Extract and count mutations
	for idx, row in df.iterrows():
		if row['Classification'] != 'WT':
			parsed = row['Mutation_info'].split(",")
			for mut in parsed:
				if "Ins" in str(mut):
					length = int(re.search(r"Ins_(\d+)", mut).group(1))
					if length not in insertion_lengths_count:
						insertion_lengths_count[length] = 0
					insertion_lengths_count[length] += 1
	
	# Create the plot
	fig = go.Figure()
	
	# Add the bars for insertion length counts
	fig.add_trace(go.Bar(
		x=list(insertion_lengths_count.keys()), 
		y=list(insertion_lengths_count.values()), 
		marker=dict(color='red'),
		width=0.8,
		opacity=1,
		name='Insertion Length Count',
		hoverinfo='x+y',
		showlegend=False
	))
	
	# Customize the layout
	fig.update_layout(
		title='Count of Insertions vs. Insertion Length',
		xaxis_title='Insertion Length',
		yaxis_title='Count of Insertions (Log Scale)',
		showlegend=True,
		yaxis_type='log',  # Set y-axis to logarithmic scale
		bargap=0.1,  # Adjust the gap between bars
		margin=dict(l=0, r=0, t=50, b=50),
		xaxis_rangeslider_visible=True	# Adjust margins if needed
	)
	
	# Show the plot
	fig.write_image(f'{output_dir}/insertion_count_length.png')
	fig.write_image(f'{output_dir}/insertion_count_length.pdf')


def allele_plot(ref_seq, cv_pos, cv_pos_2, strand, strand_2, input_file, output_dir, cleavage_pos, target, target_2, original_target, min_read_cnt, min_read_freq, plot_window, show_line, induced_mutation_str, show_all_between_allele, group_separate=False):
		
	mut_dict = []
	all_cnt = 0
	for line in open(input_file).readlines()[1:]:
		line = line.strip().split('\t')
		mut_dict.append([line[5], [int(line[1]), float(line[2]), line[3], line[4], line[0]]])
	
	mut_list = []
	for i in mut_dict:
		mut_list.append([i[0], i[1][1], i[1][0], i[1][2].split(','), i[1][3], i[1][4]])

	xlim_min = 0

	fig_wide = (plot_window * 2)* 0.185 + 15
	if cv_pos_2 != False:
		if show_all_between_allele:
			fig_wide = (plot_window * 2 + cv_pos_2 - cv_pos) * 0.2 + 6
		else:
			fig_wide = (plot_window * 4) * 0.2 + 8

	xlim_max = 0.005+(plot_window*2 + 50)*0.021
	
	window_ed = cv_pos + plot_window
	if cv_pos_2 != False:
		if show_all_between_allele:
			window_ed = cv_pos_2 + plot_window
			xlim_max = 0.005+(plot_window * 4 + cv_pos_2 - cv_pos +1)*0.021
		else:
			xlim_max = 0.005+(plot_window*5+1)*0.021
	xlim_max = fig_wide/11

	ref_plot = ref_seq[cv_pos - plot_window: window_ed]
	ref_len = len(ref_seq)

	data = {"Sequence": [], "Type": [], "Percentage": [], "LargeDel": [], "LargeIns": []}
	
	nucleotide_colors = {
		'A': 'honeydew',
		'T': 'mistyrose',
		'G': 'lightyellow',
		'C': 'aliceblue',
		'-': 'whitesmoke',
		'N': 'whitesmoke',
		' ': 'white',
	}
	
	
	print(show_line)
	fig, ax = plt.subplots(figsize=(fig_wide, 3 * (show_line+5)/10))
	
	fw = open(f'{output_dir}/Allel_table.txt', 'w') 
	fw.write('ref_seq\tmut_seq\tRaw_count\t%\tCIGAR\tMutation_info\n')

	#for info_n, info in enumerate(sorted(output_res.items(), key=lambda x: x[1], reverse=True)[:show_line]):

	#Draw Reference
	x_pos = 0.01
	add_x = 0
	seq = ref_plot
	cleavage_pos += 1
	if cv_pos_2 and not show_all_between_allele:
		add_x = 0.021 * (plot_window * 2 + 2)	

	for nucleotide in seq:
		ax.text(x_pos, 0.9, nucleotide, fontsize=14, fontfamily='monospace',
				ha='left', va='center', bbox=dict(facecolor=nucleotide_colors[nucleotide], edgecolor='none', pad=2))
		x_pos += 0.021
	if strand == 1:
		rect = patches.Rectangle((0.005+(plot_window+1-cleavage_pos)*0.021, 0.72), len(target)*0.021, 0.08, linewidth=2, edgecolor='none', facecolor='silver', zorder=10)
		ax.add_patch(rect)
		ax.text(0.005+(plot_window+1-cleavage_pos - 6)*0.021, 0.77, f'sgRNA {original_target}', fontsize=14, ha='left', va='center', zorder=11)
	else:
		rect = patches.Rectangle((0.005+(plot_window+1-(len(target) - cleavage_pos))*0.021, 0.72), len(target)*0.021, 0.08, linewidth=2, edgecolor='none', facecolor='silver', zorder=10)
		ax.add_patch(rect)
		ax.text(0.005+(plot_window+1-(len(target) - cleavage_pos) - 6)*0.021, 0.77, f'sgRNA {3 - original_target}', fontsize=14, ha='left', va='center', zorder=11)

	status_list = [[],[],[]]

	if cv_pos_2:

		if not show_all_between_allele:
			ref_plot_2 = ref_seq[cv_pos_2 - plot_window: cv_pos_2 + plot_window]
			add_x = 0.021 * (plot_window * 2 + 10)
			x_pos = 0.01 + add_x

			for nucleotide in ref_plot_2:
				ax.text(x_pos, 0.9, nucleotide, fontsize=14, fontfamily='monospace',
						ha='left', va='center', bbox=dict(facecolor=nucleotide_colors[nucleotide], edgecolor='none', pad=2))
				x_pos += 0.021 # Adjust spacing between nucleotides
			
			ax.text(0.021 * (plot_window * 2 + 5.5), 0.9, f'{cv_pos_2 - cv_pos - plot_window}bp', fontsize=10, backgroundcolor='white', ha='center', va='center', zorder=11,bbox=dict(facecolor='white', edgecolor='none', pad=2))
			target2_rec_cp = plot_window + 1

		else:
			add_x = 0
			target2_rec_cp = plot_window + cv_pos_2 - cv_pos + 1

		if strand_2 == 1:
			rect = patches.Rectangle((0.005+add_x+(target2_rec_cp-cleavage_pos)*0.021, 0.728), len(target_2)*0.021, 0.08, linewidth=2, edgecolor='none', facecolor='silver', zorder=10)
			ax.add_patch(rect)
			ax.text(0.005+add_x+(target2_rec_cp-cleavage_pos-6)*0.021, 0.77, f'sgRNA 2', fontsize=14, ha='left', va='center', zorder=11)
		else:
			rect = patches.Rectangle((0.005+add_x+(target2_rec_cp-(len(target_2) - cleavage_pos))*0.021, 0.72), len(target)*0.021, 0.08, linewidth=2, edgecolor='none', facecolor='silver', zorder=10)
			ax.add_patch(rect)
			ax.text(0.005+add_x+(target2_rec_cp-(len(target_2) - cleavage_pos) - 6)*0.021, 0.77, f'sgRNA 2', fontsize=14, ha='left', va='center', zorder=11)

	if show_all_between_allele:
		window_ed = plot_window * 2 + (cv_pos_2 - cv_pos)
	else:
		window_ed = plot_window * 2
		
	ax.text(0.021 * window_ed + 0.02 + add_x, 0.9, f'Reference', fontsize=14, ha='left', va='center')
	if induced_mutation_str:
		ax.text(0.021 * window_ed + 0.02 + add_x, 0.6, f'Desired', fontsize=14, ha='left', va='center')
	showing_line = mut_list[:show_line]

	if induced_mutation_str:
		showing_line = [[induced_mutation_str, '', '', [''], '', '']] + showing_line

	y_pos = 0.9 - 2 * 0.1
	group_separate_interval = 0

	for info_n, info in enumerate(showing_line):

		seq = ref_plot
		if group_separate and info_n != 0 and info[5] != showing_line[info_n - 1][5]:
			group_separate_interval += 1
			y_pos -= 0.1
			status_list[0].append(f'')
			status_list[1].append('')
			status_list[2].append('')
		
		x_pos = 0.01
		y_pos -= 0.1
		
		adjust_ins = 0

		if info[0] == 'None':
			for nucleotide in seq:
				ax.text(x_pos, y_pos, nucleotide, fontsize=14, fontfamily='monospace',
					ha='left', va='center', bbox=dict(facecolor=nucleotide_colors[nucleotide], edgecolor='none', pad=2))
				x_pos += 0.021 # Adjust spacing between nucleotides
			fw.write(f'{ref_plot}\t{ref_plot}\t{info[2]}\t{info[1]}\t{plot_window*2}M\t{info[0]}\n')
			status_list[0].append(f'{info[1]}% ({info[2]} Reads)')
			status_list[1].append('WT')
			status_list[2].append(info[5])
			continue
		
		mut_str = ''
	
		align_ref = ''
		align_mut = ''
		CIGAR_str = ''
		
		info_sp = info[0].split(',')

		ex_mut = info_sp[0].replace(':', '_').replace('>', '_').split('_')
		ex_mut[1] = int(ex_mut[0]) - plot_window - 1

		right_large_check = False
		right_large_ins_check = False
		left_large_check = False
		left_large_ins_check = False

		for mut_n, mut in enumerate(info_sp):
			
			mut = mut.replace(':', '_').replace('>', '_').split('_')
			for i_n, i in enumerate(mut):
				if i.isdigit():
					mut[i_n] = int(i)

			align_ref += ref_seq[ex_mut[1] + 1: mut[0]]
			align_mut += ref_seq[ex_mut[1] + 1: mut[0]]

			dist = mut[0] - ex_mut[1] - 1
			if dist > 0:
				CIGAR_str += str(mut[0] - ex_mut[1] - 1) + 'M'

			if mut[2] == 'Sub':
				align_ref += mut[4]
				align_mut += mut[5]
				CIGAR_str += f'{mut[3]}='
			elif 'Del' in mut[2]:
				align_ref += ref_seq[mut[0]: mut[1] + 1]
				align_mut += '-'*mut[3]
				CIGAR_str += f'{mut[3]}D'
				if int(mut[0]) - (cv_pos - plot_window) + adjust_ins < 0:
					left_large_check = True
					left_seq = ref_seq[mut[0] - 10: mut[0]]
					#left_len_comp = 2 * 0.021 * (1- (mut[0]/ cv_pos_2))
					left_len_comp = 0
					xlim_min = -0.005-(11+2)*0.021
				if int(mut[1]) - (cv_pos - plot_window) + adjust_ins > window_ed:
					right_large_check = True
					right_seq = ref_seq[mut[1]-10: mut[1]]
					right_len_comp = 0
				if (right_large_check or left_large_check) and len(info_sp) > mut_n + 1:
						next_mut = info_sp[mut_n + 1].replace(':', '_').replace('>', '_').split('_')
						if mut[0] == int(next_mut[0]) and 'Ins' in next_mut[2]:
							ins_len = int(next_mut[3])
							if right_large_check:
								right_large_ins_check = True
							if left_large_check:
								left_large_ins_check = True
			elif 'Ins' in mut[2]:
				align_ref += '-'*mut[3]
				align_mut += mut[4]
				CIGAR_str += f'{mut[3]}I'

			ex_mut = mut
		
		align_ref += ref_seq[mut[1] + 1: mut[1] + plot_window]
		align_mut += ref_seq[mut[1] + 1: mut[1] + plot_window]
		CIGAR_str += f'{plot_window}M'

		if not induced_mutation_str or info_n != 0:
			fw.write(f'{align_ref}\t{align_mut}\t{info[2]}\t{info[1]}\t{CIGAR_str}\t{info[0]}\n')

		bar_st_between = 0
		bar_ed_between = 10

		inv_info = False
		for i in info[3]:
			if i[:9] == 'inversion':
				i = i[i.find('_')+1:i.find('inSeq')].split('to')
				inv_info = [int(i[0]), int(i[1])]

		for mut_n, mut in enumerate(info_sp):
			
			mut_st = int(mut.split(':')[0].split('_')[0])
			mut_ed = int(mut.split(':')[0].split('_')[1])

			if mut.find('Sub') != -1:
				st_pos = mut_st - (cv_pos - plot_window) + adjust_ins
				ed_pos = len(ref_seq)
			else:
				st_pos = mut_st - (cv_pos - plot_window) + adjust_ins
				ed_pos = mut_ed - (cv_pos - plot_window) + adjust_ins
		
	
			mut_type = mut.split(':')[1].split('_')[0]
			
			if st_pos > window_ed:
				if mut_type == 'Sub':
					mut_str += f'{mut.split(":")[1].split("_")[1]}Sub,'
				else:
					mut_len = int(mut.split(':')[1].replace('Large','Large').split('_')[1])
					if mut_type == 'Del':
						mut_str += f'{mut_len}Del,'
					elif mut_type == 'Ins':
						mut_str += f'{mut_len}Ins,'
					elif mut_type == 'LargeDel':
						mut_str += f'{mut_len}LargeDel,'
						if not show_all_between_allele:
							if cv_pos_2:
								add_len_between = (mut_ed - cv_pos - plot_window) / (cv_pos_2 - cv_pos) * 10
								if add_len_between > 10:
									add_len_between = 10
								if bar_st_between < add_len_between:
									bar_st_between = add_len_between
								#ax.annotate('', xy=(0.005+60*0.021, y_pos), xytext=(0.01+(ed_pos+1+add_len_between)*0.021, y_pos), arrowprops=dict(arrowstyle='-'), zorder=10)
							if mut_ed > cv_pos_2 - plot_window and mut_st < cv_pos_2 - plot_window:
								add_len_between = (mut_st - cv_pos - plot_window) / (cv_pos_2 - cv_pos) * 10
								if add_len_between > 10:
									add_len_between = 10
								if bar_ed_between > add_len_between:
									bar_ed_between = add_len_between

					elif mut_type == 'LargeIns':
						mut_str += f'{mut_len}LargeIns,'

				continue
			
			out_of_range = False
			if st_pos < 0:
				st_pos = 0
			if ed_pos >= window_ed:
				ed_pos = window_ed - 1
				out_of_range = True
			
			if mut_type == 'Del':
				del_len = int(mut.split(':')[1].split('_')[1])
				if out_of_range == False:
					seq = seq[:st_pos] + '-'*(ed_pos - st_pos + 1) + seq[ed_pos+1:]
				mut_str += f'{del_len}Del,'
			
			elif mut_type == 'Ins':
				ins_len = int(mut.split(':')[1].split('_')[1])
				if out_of_range == False:
					rect = patches.Rectangle((0.008+st_pos*0.021,	y_pos -0.04), ins_len*0.021, 0.1, linewidth=2, edgecolor='r', facecolor='none', zorder=10)
					ax.add_patch(rect)
					seq = (seq[:st_pos] + mut.split('_')[3] + seq[st_pos:])[:window_ed]
					adjust_ins += ins_len
				mut_str += f'{ins_len}Ins,'
			
			elif mut_type == 'LargeDel':
				del_len = int(mut.split(':')[1].replace('Large','Large').split('_')[1])
				seq = seq[:st_pos] + ' '*(ed_pos - st_pos + 1) + seq[ed_pos+1:]
				inv_str = ''
				if inv_info:
					if inv_info[0] - 10 < mut_st  and mut_ed < inv_info[1] + 10:
						inv_str = ' (inv)'
				if cv_pos_2 and not show_all_between_allele:
					add_len_between = (mut_ed - cv_pos - plot_window) / (cv_pos_2 - cv_pos  - plot_window*2) * 10
					if add_len_between > 10:
						add_len_between = 10
					if bar_st_between < add_len_between:
						bar_st_between = add_len_between
					add_large_cv2 = 10.5
				else:
					add_len_between = 0
					add_large_cv2 = 0
				if left_large_check or right_large_check:
					if left_large_check and right_large_check:
						ax.annotate('', xy=(-0.005-left_len_comp, y_pos), xytext=(0.005+(window_ed+add_large_cv2)*0.021+right_len_comp, y_pos), arrowprops=dict(arrowstyle='-'), zorder=9)
					elif left_large_check:
						ax.annotate('', xy=(-0.005-left_len_comp, y_pos), xytext=(0.01+(ed_pos+1+add_len_between+0.05)*0.021, y_pos), arrowprops=dict(arrowstyle='-'), zorder=9) 	
					elif right_large_check:
						ax.annotate('', xy=(0.005+(st_pos+adjust_ins)*0.021, y_pos), xytext=(0.005+(window_ed+add_large_cv2)*0.021+right_len_comp, y_pos), arrowprops=dict(arrowstyle='-'), zorder=9) 	
					if left_large_ins_check:
						rect = patches.Rectangle((-0.005-left_len_comp, y_pos -0.04), 0.005, 0.1, linewidth=2, edgecolor='r', facecolor='none', zorder=9)
						ax.add_patch(rect)
					if right_large_ins_check and cv_pos_2 == False:
						rect = patches.Rectangle((0.005+window_ed*0.021+right_len_comp,	y_pos -0.04), 0.006, 0.1, linewidth=2, edgecolor='r', facecolor='none', zorder=13)
						ax.add_patch(rect)
						mut_str += f'{ins_len}LargeIns,'
				else:
					ax.annotate('', xy=(0.005+(st_pos+adjust_ins)*0.021, y_pos), xytext=(0.05+(ed_pos+1+add_len_between+0.01)*0.021, y_pos), arrowprops=dict(arrowstyle='-'), zorder=10)
				if cv_pos_2 and not show_all_between_allele:
					if mut_ed > cv_pos_2 - plot_window:
						ax.text(0.005+(plot_window*2+5)*0.021, y_pos, f'{del_len}bp{inv_str}', fontsize=10, backgroundcolor='white', ha='center', va='center', zorder=11,bbox=dict(facecolor='white', edgecolor='none', pad=2))		
					else:
						ax.text(0.005+(st_pos+ed_pos + add_len_between)/2*0.021, y_pos, f'{del_len}bp{inv_str}', fontsize=10, backgroundcolor='white', ha='left', va='center', zorder=11,bbox=dict(facecolor='white', edgecolor='none', pad=2))	
				else:
					ax.text(0.005+(st_pos+ed_pos)/2*0.021, y_pos, f'{del_len}bp{inv_str}', fontsize=10, backgroundcolor='white', ha='left', va='center', zorder=11,bbox=dict(facecolor='white', edgecolor='none', pad=2))
				mut_str += f'{del_len}LargeDel,'
				
			elif mut_type == 'LargeIns':
				ins_len = int(mut.split(':')[1].replace('Large','Large').split('_')[1])
				mut_str += f'{ins_len}LargeIns,'
				if right_large_ins_check or left_large_ins_check:
					if int(info_sp[mut_n - 1].split(':')[0].split('_')[0]) == mut_st:
						continue
				#seq = (seq[:st_pos] + mut.split('_')[3] + seq[st_pos:])[:plot_window*2]
				if st_pos + ins_len > window_ed:
					ins_len = window_ed - st_pos
				rect = patches.Rectangle((0.008+st_pos*0.021,	y_pos -0.04), ins_len*0.021, 0.1, linewidth=2, edgecolor='r', facecolor='none', zorder=14) 
				ax.add_patch(rect)
				adjust_ins += ins_len
			
			elif mut_type == 'Sub':
				mut_nt = mut.split('>')[1]
				seq = seq[:st_pos] + mut_nt.lower() + seq[st_pos+len(mut_nt):]
				mut_str += f'{mut.split(":")[1].split("_")[1]}Sub,'
			
		if info[3] != ['']:
			for insert in info[3][0].split(','):
				mut_str += f'{insert.split("_")[0]},'
			mut_str = mut_str
			
		seq = seq[:window_ed]

		for nucleotide in seq:
			if nucleotide in 'atgcn':
				nucleotide = nucleotide.upper()
				ax.text(x_pos, y_pos, nucleotide, fontsize=14, fontfamily='monospace', fontweight='bold',
					ha='left', va='center', bbox=dict(facecolor=nucleotide_colors[nucleotide], edgecolor='none', pad=2))
			else:
				ax.text(x_pos, y_pos, nucleotide, fontsize=14, fontfamily='monospace',
					ha='left', va='center', bbox=dict(facecolor=nucleotide_colors[nucleotide], edgecolor='none', pad=2))
			x_pos += 0.021 # Adjust spacing between nucleotides

		if left_large_check:
			left_x_pos = -0.021 * 1 - left_len_comp
			for nucleotide in left_seq[::-1]:
				ax.text(left_x_pos, y_pos, nucleotide, fontsize=14, fontfamily='monospace',
					ha='left', va='center', bbox=dict(facecolor=nucleotide_colors[nucleotide], edgecolor='none', pad=2))
				left_x_pos -= 0.021 # Adjust spacing between nucleotides
		
		if right_large_check and (not cv_pos_2 or (cv_pos_2 and show_all_between_allele)):
			right_x_pos = 0.01 + 0.021 * window_ed + right_len_comp
			for nucleotide in right_seq:
				ax.text(right_x_pos, y_pos, nucleotide, fontsize=14, fontfamily='monospace',
					ha='left', va='center', bbox=dict(facecolor=nucleotide_colors[nucleotide], edgecolor='none', pad=2))
				right_x_pos += 0.021 # Adjust spacing between nucleotides
		
		if induced_mutation_str and info_n == 0:
			status_list[0].append(f'')
			status_list[1].append('')
			#status_list[2].append('')
			status_list[2].append(info[5])
		else:
			status_list[0].append(f'{info[1]}% ({info[2]} Reads)')
			status_list[1].append(info[4])
			#status_list[2].append(mut_str[:-1])
			status_list[2].append(info[5])
			
		if cv_pos_2 and bar_st_between < bar_ed_between and not show_all_between_allele:
			rect = patches.Rectangle((0.0075+ plot_window*2*0.021 + bar_st_between,y_pos -0.035), (bar_ed_between - bar_st_between)*0.021, 0.08, linewidth=0, edgecolor='none', facecolor='whitesmoke', zorder=13)
			ax.add_patch(rect)

		if induced_mutation_str and info_n == 0:
			y_pos -= 0.15

	if len(mut_list) < show_line:
		dash_line = len(mut_list)
	else:
		dash_line = show_line

	if induced_mutation_str:
		dash_line += 2
	dash_line += group_separate_interval

	st_pos = 0
	for i in range(len(status_list)):
		max_len = 0
		for info_n, x in enumerate(status_list[i]):
			y_pos = 0.9 - (info_n + 3) * 0.1
			if induced_mutation_str:
				y_pos -= 0.15
			ax.text(0.021 * window_ed + 0.02 + add_x + (11+2 +st_pos)*0.021, y_pos, x, fontsize=12, ha='left', va='center')
			if len(x) > max_len:
				max_len = len(x)
		st_pos += max_len/2 + 2

	ax.annotate('', xy=(0.005+(plot_window)*0.021, 1), xytext=(0.01+(plot_window)*0.021, 1 - dash_line/10 -0.4), arrowprops=dict(arrowstyle='-', linestyle='--', color='gray'), zorder=10)
	ax.text(0.021 * plot_window * 2 + 1, y_pos, f'  ', fontsize=12, ha='left', va='center')

	if cv_pos_2 and not show_all_between_allele:
		
		y_pos = 0.9 - 2 * 0.1

		for info_n, info in enumerate(showing_line):

			seq = ref_plot_2
			
			x_pos = 0.01 + add_x
			y_pos -= 0.1
			
			if group_separate and info_n != 0 and info[5] != showing_line[info_n - 1][5]:
				y_pos -= 0.1
			
			
			adjust_ins = 0
			
			if info[0] == 'None':
				for nucleotide in seq:
					ax.text(x_pos, y_pos, nucleotide, fontsize=14, fontfamily='monospace',
							ha='left', va='center', bbox=dict(facecolor=nucleotide_colors[nucleotide], edgecolor='none', pad=2))
					x_pos += 0.021 # Adjust spacing between nucleotides
				fw.write(f'{ref_plot_2}\t{ref_plot_2}\t{info[2]}\t{info[1]}\t{plot_window*2}M\t{info[0]}\n')
				rect = patches.Rectangle((0.0075+plot_window*2*0.021,y_pos -0.035), 10 * 0.021, 0.08, linewidth=0, edgecolor='none', facecolor='whitesmoke', zorder=13)
				ax.add_patch(rect)
				continue
		
			align_ref = ''
			align_mut = ''
			CIGAR_str = ''
			
			info_sp = info[0].split(',')
			ex_mut = info_sp[0].replace(':', '_').replace('>', '_').split('_')
			ex_mut[1] = int(ex_mut[0]) - plot_window - 1

			right_large_check = False
			right_large_ins_check = False
			left_large_check = False
			left_not_draw_between = False
			left_large_ins_check = False

			for mut_n, mut in enumerate(info_sp):
				mut = mut.replace(':', '_').replace('>', '_').split('_')

				for i_n, i in enumerate(mut):
					if i.isdigit():
						mut[i_n] = int(i)

				align_ref += ref_seq[ex_mut[1] + 1: mut[0]]
				align_mut += ref_seq[ex_mut[1] + 1: mut[0]]

				dist = mut[0] - ex_mut[1] - 1
				if dist > 0:
					CIGAR_str += str(mut[0] - ex_mut[1] - 1) + 'M'

				if mut[2] == 'Sub':
					align_ref += mut[4]
					align_mut += mut[5]
					CIGAR_str += f'{mut[3]}='
				elif 'Del' in mut[2]:
					align_ref += ref_seq[mut[0]: mut[1] + 1]
					align_mut += '-'*mut[3]
					CIGAR_str += f'{mut[3]}D'
					if int(mut[0]) - (cv_pos_2 - plot_window) + adjust_ins < 0 and int(mut[3]) > 100:
						left_large_check = True
						if mut[0] > cv_pos + plot_window:
							left_not_draw_between = True
						left_seq = ref_seq[mut[0] - 10: mut[0]]
						#right_len_comp = 2 * 0.021 * (1- (mut[0]/ cv_pos_2))
						left_len_comp = 0
					if int(mut[1]) - (cv_pos_2 - plot_window) + adjust_ins > plot_window*2:
						right_large_check = True
						right_seq = ref_seq[mut[1]: mut[1]+10]
						#left_len_comp = 2 * 0.021 * (1- ((ref_len - mut[1])/(ref_len - cv_pos_2)))
						right_len_comp = 0
					if (left_large_check or right_large_check) and len(info_sp) > mut_n + 1:
						next_mut = info_sp[mut_n + 1].replace(':', '_').replace('>', '_').split('_')
						if mut[0] == int(next_mut[0]) and 'Ins' in next_mut[2]:
							if right_large_check:
								right_large_ins_check = True
							if left_large_check:
								left_large_ins_check = True
				elif 'Ins' in mut[2]:
					align_ref += '-'*mut[3]
					align_mut += mut[4]
					CIGAR_str += f'{mut[3]}I'

				ex_mut = mut

			align_ref += ref_seq[mut[1] + 1: mut[1] + plot_window]
			align_mut += ref_seq[mut[1] + 1: mut[1] + plot_window]
			CIGAR_str += f'{plot_window}M'

			if not induced_mutation_str or info_n != 0:
				fw.write(f'{align_ref}\t{align_mut}\t{info[2]}\t{info[1]}\t{CIGAR_str}\t{info[0]}\n')


			for mut in info[0].split(','):

				mut_st = int(mut.split(':')[0].split('_')[0])
				mut_ed = int(mut.split(':')[0].split('_')[1])

				if mut.find('Sub') != -1:
					st_pos = mut_st - (cv_pos_2 - plot_window) + adjust_ins
					ed_pos = len(ref_seq)
				else:
					st_pos = mut_st - (cv_pos_2 - plot_window) + adjust_ins
					ed_pos = mut_ed - (cv_pos_2 - plot_window) + adjust_ins
		
				mut_type = mut.split(':')[1].split('_')[0]
				
				if st_pos > 60:
					continue
				
				out_of_range = False
				if st_pos < 0:
					st_pos = 0
					out_of_range = True
				if ed_pos >= plot_window*2:
					ed_pos = plot_window*2 - 1
				
				if mut_type == 'Del':
					if out_of_range == False:
						del_len = int(mut.split(':')[1].split('_')[1])
						seq = seq[:st_pos] + '-'*(ed_pos - st_pos + 1) + seq[ed_pos+1:]
					if cv_pos + plot_window <= mut_st <= cv_pos_2 - plot_window or cv_pos + plot_window <= mut_ed <= cv_pos_2 - plot_window:
						mark_st = mut_st-cv_pos-plot_window
						mark_ed = mut_ed-cv_pos-plot_window
						mark_len = cv_pos_2-cv_pos-plot_window*2
						if mark_st < 0:
							mark_st = 0
						if mark_ed > mark_len:
							mark_ed = mark_len
						ax.annotate('', xy=((plot_window*2 + mark_st*10/mark_len)*0.021, y_pos), xytext=((plot_window*2 + mark_ed*10/mark_len)*0.021, y_pos), arrowprops=dict(arrowstyle='-', color='black'), zorder=15)
				
				elif mut_type == 'Ins':
					ins_len = int(mut.split(':')[1].split('_')[1])
					if out_of_range == False:
						rect = patches.Rectangle((0.008+add_x+st_pos*0.021, y_pos - 0.04), ins_len*0.021, 0.1, linewidth=2, edgecolor='r', facecolor='none', zorder=11)
						ax.add_patch(rect)
						seq = (seq[:st_pos] + mut.split('_')[3] + seq[st_pos:])
						adjust_ins += ins_len
					if cv_pos + plot_window <= mut_st and mut_ed < cv_pos_2 - plot_window:
						ax.annotate('', xy=((plot_window*2 + (mut_st-cv_pos-plot_window)*10/(cv_pos_2-cv_pos-plot_window*2))*0.021, y_pos-0.05), xytext=((plot_window*2 + (mut_st-cv_pos-plot_window)*10/(cv_pos_2-cv_pos-plot_window*2))*0.021, y_pos+0.05), arrowprops=dict(arrowstyle='-', color='red'), zorder=15)
				
				elif mut_type == 'LargeDel' and ed_pos > 0:
					del_len = int(mut.split(':')[1].replace('Large','Large').split('_')[1])
					seq = seq[:st_pos] + ' '*(ed_pos - st_pos + 1) + seq[ed_pos+1:]
					add_len_between = 0
					if inv_info:
						if inv_info[0] - 10 < mut_st  and mut_ed  < inv_info[1] + 10:
							inv_str = ' (inv)'					
					if  cv_pos_2 - plot_window > mut_st > cv_pos + plot_window:
						add_len_between = (mut_st - cv_pos - plot_window) / (cv_pos_2 - cv_pos - plot_window*2) * 10

					if right_large_check:
						ax.annotate('', xy=(0.005+add_x+st_pos*0.021, y_pos), xytext=(0.01+add_x+(ed_pos+1)*0.021 + right_len_comp, y_pos), arrowprops=dict(arrowstyle='-'), zorder=11)
					else:
						ax.annotate('', xy=(0.005+add_x+(st_pos - add_len_between)*0.021, y_pos), xytext=(0.01+add_x+(ed_pos+1)*0.021, y_pos), arrowprops=dict(arrowstyle='-'), zorder=11)
					if right_large_ins_check:
						rect = patches.Rectangle((0.01+add_x+(ed_pos+1)*0.021 + right_len_comp- 0.006, y_pos -0.04), 0.006, 0.1, linewidth=2, edgecolor='r', facecolor='none', zorder=13)
						ax.add_patch(rect)
					if cv_pos_2 - plot_window < mut_st:
						ax.text(0.005+add_x+(st_pos-add_len_between+ed_pos)/2*0.021, y_pos, f'{del_len}bp{inv_str}', fontsize=10, backgroundcolor='white', ha='left', va='center', zorder=11,bbox=dict(facecolor='white', edgecolor='none', pad=2))
					if cv_pos + plot_window <= mut_st and mut_ed < cv_pos_2 - plot_window:
						ax.annotate('', xy=((plot_window*2 + (mut_st-cv_pos-plot_window)*10/(cv_pos_2-cv_pos-plot_window*2))*0.021, y_pos), xytext=((plot_window*2 + (mut_ed-cv_pos-plot_window)*10/(cv_pos_2-cv_pos-plot_window*2))*0.021, y_pos), arrowprops=dict(arrowstyle='-', color='black'), zorder=15)
					
				elif mut_type == 'LargeIns':
					if right_large_ins_check or left_large_ins_check:
						if int(info_sp[mut_n - 1].split(':')[0].split('_')[0]) == mut_st:
							continue
					ins_len = int(mut.split(':')[1].replace('Large','Large').split('_')[1])
					seq = (seq[:st_pos] + mut.split('_')[3] + seq[st_pos:])[:plot_window*2]
					if st_pos + ins_len > plot_window *2:
						ins_len = plot_window*2 - st_pos
					rect = patches.Rectangle((0.008+add_x+st_pos*0.021, y_pos - 0.04), ins_len*0.021, 0.1, linewidth=2, edgecolor='r', facecolor='none', zorder=13)
					ax.add_patch(rect)
					adjust_ins += ins_len
					if cv_pos + plot_window <= mut_st and mut_ed < cv_pos_2 - plot_window:
						ax.annotate('', xy=((plot_window*2 + (mut_st-cv_pos-plot_window)*10/(cv_pos_2-cv_pos-plot_window*2))*0.021, y_pos-0.05), xytext=((plot_window*2 + (mut_st-cv_pos-plot_window)*10/(cv_pos_2-cv_pos-plot_window*2))*0.021, y_pos+0.05), arrowprops=dict(arrowstyle='-', color='red'), zorder=15)
				
				elif mut_type == 'Sub':
					for n, nt in enumerate(mut.split('>')[1]):
						if cv_pos_2 - plot_window <= mut_st + n < cv_pos_2 + plot_window:
							seq = seq[:mut_st - (cv_pos_2 - plot_window)] + nt.lower() + seq[mut_st - (cv_pos_2 - plot_window)+1:]
						if cv_pos + plot_window <= mut_st + n < cv_pos_2 - plot_window:
							ax.annotate('', xy=((plot_window*2 + (mut_st-cv_pos-plot_window)*10/(cv_pos_2-cv_pos-plot_window*2))*0.021, y_pos-0.05), xytext=((plot_window*2 + (mut_st-cv_pos-plot_window)*10/(cv_pos_2-cv_pos-plot_window*2))*0.021, y_pos+0.05), arrowprops=dict(arrowstyle='-', color='black'), zorder=15)		

			seq = seq[:plot_window*2]
					
			for nucleotide in seq:
				if nucleotide in 'atgcn':
					nucleotide = nucleotide.upper()
					ax.text(x_pos, y_pos, nucleotide, fontsize=14, fontfamily='monospace', fontweight='bold',
						ha='left', va='center', bbox=dict(facecolor=nucleotide_colors[nucleotide], edgecolor='none', pad=2))
				else:
					ax.text(x_pos, y_pos, nucleotide, fontsize=14, fontfamily='monospace',
						ha='left', va='center', bbox=dict(facecolor=nucleotide_colors[nucleotide], edgecolor='none', pad=2))
				x_pos += 0.021 # Adjust spacing between nucleotides

			if right_large_check:
				x_pos += right_len_comp
				for nucleotide in right_seq:
					ax.text(x_pos, y_pos, nucleotide, fontsize=14, fontfamily='monospace',
							ha='left', va='center', bbox=dict(facecolor=nucleotide_colors[nucleotide], edgecolor='none', pad=2))
					x_pos += 0.021

			if left_large_check and left_not_draw_between:
				rect = patches.Rectangle((0.0075+plot_window*2*0.021,y_pos -0.035), (10 - add_len_between)*0.021, 0.08, linewidth=0, edgecolor='none', facecolor='whitesmoke', zorder=13)
				ax.add_patch(rect)

			if induced_mutation_str and info_n == 0:
				y_pos -= 0.15

	if cv_pos_2:
		if show_all_between_allele:
			cv_pos_2_plot = 0.01+(plot_window + cv_pos_2 - cv_pos + 1)*0.021
		else:
			cv_pos_2_plot = 0.01+add_x+(plot_window+1)*0.021
		ax.annotate('', xy=(cv_pos_2_plot, 1), xytext=(cv_pos_2_plot, 1 - dash_line/10-0.4), arrowprops=dict(arrowstyle='-', linestyle='--', color='red'), zorder=10)

	ax.text(0.021 * plot_window * 2 + 1 + add_x, y_pos, f'  ', fontsize=12, ha='left', va='center')

	plt.xlim(xlim_min,xlim_max)
	plt.ylim(y_pos-0.1, 1 )
	plt.axis('off')
	plt.savefig(output_dir + '/allel_plot.png', bbox_inches='tight')
	plt.savefig(output_dir + '/allel_plot.pdf', bbox_inches='tight')


def LD_tornado(tsv_file, cv_pos, ref_len, strand, output_dir):
		
	fig, ax = plt.subplots(figsize=(12, 10))
	
	LD_list =[]

	for line in open(tsv_file).readlines():
		line = line.split('\t')
		if line[5].find('LargeDel') == -1:
			continue
		LD_st = int(line[5].split('_')[0])
		LD_ed = int(line[5].split(':')[0].split('_')[1])
		for i in range(int(line[1])):
			LD_list.append([LD_st, LD_ed, LD_ed - LD_st + 1])

	if len(LD_list) > 20:
		y_interval = 2.5 / len(LD_list)
	else:
		y_interval = 0.1
	
	plt.xlim(0,ref_len)
	plt.ylim(-1, 2)
	# Add annotations for deletions and insertions
	
	if strand == 1:
		rect = patches.Rectangle((cv_pos-ref_len*0.05, 1.92), ref_len*0.06, 0.05, linewidth=1, edgecolor='black', facecolor='lightgray', zorder=10)
		ax.add_patch(rect)
		ax.text(cv_pos+ref_len*0.014 , 1.92, 'PAM', fontsize=12, backgroundcolor='white', ha='left', zorder=9)
	elif strand == -1:
		rect = patches.Rectangle((cv_pos-ref_len*0.01, 1.92), ref_len*0.06, 0.05, linewidth=1, edgecolor='black', facecolor='lightgray', zorder=10)
		ax.add_patch(rect)
		ax.text(cv_pos-ref_len*0.014 , 1.92, 'PAM', fontsize=12, backgroundcolor='white', ha='right', zorder=9)
	ax.annotate('', xy=(0, 1.945), xytext=(ref_len, 1.945), arrowprops=dict(arrowstyle='-'))
	
	y_val = 1.75
	for mut_info in sorted(LD_list, key=lambda x: x[2], reverse=True):
	
		y_val -= y_interval
		if mut_info[0] < cv_pos-100 and mut_info[1] > cv_pos+100:
			line_color = 'gray'
		elif abs(mut_info[0]-cv_pos) < abs(mut_info[1]-cv_pos):
			line_color = 'blue'
		elif abs(mut_info[0]-cv_pos) > abs(mut_info[1]-cv_pos):
			line_color = 'red'
		
		ax.annotate('', xy=(mut_info[0], y_val), xytext=(mut_info[1], y_val), arrowprops=dict(arrowstyle='-', linewidth=1, color=line_color))
		
	ax.annotate('', xy=(cv_pos, 1.8), xytext=(cv_pos, -1), arrowprops=dict(arrowstyle='-', linestyle='--', linewidth=1, color='gray'))
	
	# Add a vertical dotted line
	
	ax.xaxis.set_ticks_position('top')
	ax.xaxis.set_label_position('top')
	ax.yaxis.set_visible(False)
	
	ax.spines['top'].set_position(('outward', -40))
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.spines['right'].set_visible(False)
	
	plt.savefig(output_dir + '/Largedeletion_tornado.png')
	plt.savefig(output_dir + '/Largedeletion_tornado.pdf')

def write_html(plots, control_check, target_check, output_dir):

	full_html = """
	<!DOCTYPE html>
	<html>
	<head>
		<title>Combined Graphs</title>
		<style>
			body {margin: 0;padding: 0;display: flex;justify-content: center;align-items: center;background-color: #f4f4f4;}
    		#graph-container {min-width: 500px; max-width: 1000px;width: 80vw;background-color: #ffffff;border: 1px solid #ddd;box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);overflow: hidden; }
			img {max-width: 100%;max-height: 100%;object-fit: contain;}
		</style>
	</head>
	<body>
		<div id='graph-container'>"""
	
	full_html += f"<h1>Preprocess graphs</h1>{plots['treated_align']}{plots['control_align']}<br />"

	if control_check:
		full_html += f'<h1>Accuracy graphs</h1><img src="results/Regular_accuracy.png" alt="Matplotlib Graph"/><br />'

	full_html += f'<h1>Allele plot</h1><img src="results/allel_plot.png" alt="Matplotlib Graph" /><br />'
	full_html += f'<h1>Mutation plot</h1>{plots["mutation_pie"]}{plots["pattern_pie"]}{plots["allele_pie"]}<br />'
	full_html += f'<h1>Indel plot</h1>{plots["indel_per_pos"]}<br />'
	full_html += f'<h1>Insertion plot</h1>{plots["insertion_len"]}<br />'
	full_html += f'<h1>Deletion plot</h1>{plots["deletion_len"]}{plots["deletion_count_len"]}<br />'
	full_html += f'<h1>Large Deletion plot</h1><img src="results/Largedeletion_tornado.png" alt="Matplotlib Graph" /><br />'
	full_html += f'<h1>Substitutions plot</h1>{plots["base_proportion"]}<br />'


	full_html += """</div></body>
	</html>
	"""

	with open(output_dir + "/combined_graphs.html", "w") as fw:
		fw.write(full_html)


def write_read_count(tsv_file, input_pre_cnt_file, output_read_file, output_summary_file, min_read_cnt, min_read_freq, induced_sequence_path, custom_all_cnt = False):

	mut_dict = {}
	all_cnt = 0
	for line in open(tsv_file).readlines()[1:]:
		line_sp = line.strip().split('\t')
		mut = line_sp[2]
		all_cnt += 1
		if mut not in mut_dict:
			if len(line_sp) < 5:
				mut_dict[mut] = [1, '']
			else: 
				mut_dict[mut] = [1, line_sp[3], line_sp[4], line_sp[1]]
		else:
			mut_dict[mut][0] += 1
	
	filt_mut_dict= {}
	filt_all_cnt = 0
	filted_low_cnt = 0
	for i in sorted(mut_dict.items(), key=lambda x: x[1][0], reverse=True):
		if i[1][0] < min_read_cnt or round(i[1][0]*100/all_cnt,2) < min_read_freq:
			filted_low_cnt += i[1][0]
			continue
		filt_mut_dict[i[0]] = i[1]
		filt_all_cnt += i[1][0]

	mut_dict = filt_mut_dict
	all_cnt = filt_all_cnt

	mut_cnt = {}
	for i in ['WT', 'Del', 'Ins', 'Sub', 'LargeDel', 'LargeIns', 'Inv','Complex', 'ComplexWithLargeMut']:
		mut_cnt[i] = 0
	precise_cnt = {}
	for i in ['Precise', 'Precise_with_mutations', 'Partial', 'Partial_with_mutations', 'X']:
		precise_cnt[i] = 0

	with open(output_read_file, 'w') as fw:
		fw.write('Mutation_patterns\tCount\t(%)\tIntegration_information\tInduced_mutation\tDetail_mutation_informations\n')
		for i in mut_dict.items():
			fw.write(f'{i[1][3]}\t{i[1][0]}\t{round(i[1][0]*100/all_cnt,2)}\t{i[1][1][:-1]}\t{i[1][2]}\t{i[0]}\n')
			if i[1][3] not in mut_cnt:
				mut_cnt[i[1][3]] = 1
			else:
				mut_cnt[i[1][3]] += i[1][0]
			if i[1][2] not in precise_cnt:
				precise_cnt[i[1][2]] = 1
			else:
				precise_cnt[i[1][2]] += i[1][0]
			

	pre_dict = {}
	f = open(input_pre_cnt_file).readlines()
	cnt = f[1].strip().split('\t')
	for i in zip(f[0].strip().split('\t'), cnt):
		pre_dict[i[0]] = i[1]

	fw = open(output_summary_file, 'w') 
	s1 = ''
	s2 = ''

	fw.write('\n')
	print(pre_dict)
	for i in ['Treated_all_reads','Treated_unmapped','Treated_low_quality','Treated_short', 'all', 'unmapped', 'filted_align_short', 'filted_umi_len', 'used', 'used_with_supple']:
		if i  in pre_dict:
			s1 += i + '\t'
			s2 += str(pre_dict[i]) + '\t'
	if 'Treated_used' in pre_dict:
		s1 += 'Treated_low_cnt_reads\tTreated_used\t'
		s2 += str(filted_low_cnt) + '\t' + str(int(pre_dict['Treated_used']) - filted_low_cnt) + '\t'
	
	if custom_all_cnt:
		s1 += 'Custom_classified_read\t'
		s2 += int(custom_all_cnt) + '\t'

	s1 += '\t'
	s2 += '\t'

	for i in ['WT', 'Del', 'Ins', 'Sub', 'LargeDel', 'LargeIns', 'Inv', 'Complex', 'ComplexWithLargeMut']:
		if i  in mut_cnt:
			s1 += i + '\t'
			s2 += str(mut_cnt[i]) + '\t'
			print(i, mut_cnt[i])
	s1 += '\t'
	s2 += '\t'

	if induced_sequence_path != False:
		for i in ['Precise', 'Precise_with_mutations', 'Partial', 'Partial_with_mutations', 'X']:
			if i  in precise_cnt:
				s1 += i + '\t'
				s2 += str(precise_cnt[i]) + '\t'

	fw.write(s1 + '\n' + s2 + '\n')
	fw.close()


	
	
	
		
	
	
