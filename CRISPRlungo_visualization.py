#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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


# In[26]:


def visualization_preprocess(sam_file, fasta_file):

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


# Read processing information

# In[27]:


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
        plot_bgcolor='white'
    )
    
    # Show the plot
    fig.write_image(f'{output_dir}/Alignment_statistics.png')


# Create interactive graph of reads per base and show porportion of each base at each position.

# In[28]:


def base_proportion(result, output_dir, reference):

    # Assume result is already defined and populated as in your initial script
    
    # Calculate total reads and base proportions
    summed_results = {}
    base_proportions = {'A': [], 'T': [], 'G': [], 'C': []}
    
    for key in result.keys():
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
    base_colors = {'A': '#7fc97f', 'T': '#beaed4', 'G': '#ffff99', 'C': '#fdc086'}
    
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
            y=[p * maxreads/10 for  p in base_proportions[base]], 
            base=[(-maxreads/10)+adj for adj in ajustment],
            marker_color=base_colors[base], 
            marker_line_width=0, 
            name=f'{base} Proportion',
            width=0.8, 
            offset=0.1,  
        ))
        ajustment=[adj+p * maxreads/10 for adj,p in zip(ajustment,base_proportions[base])]
    
    
    for base in ['A', 'T', 'G', 'C']:
        fig.add_trace(go.Bar(
            x=np.add(ref_seq_dictionary[base], -0.5),
            y=[maxreads/10 for a in ref_seq_dictionary[base]], 
            base=[(-maxreads/4) for a in ref_seq_dictionary[base]],
            marker_color=base_colors[base], 
            showlegend=False,
            marker_line_width=0, 
            width=0.8, 
            offset=0.1,  
        ))
        
    fig.update_layout(
        title='Reads by Position with Base Proportions',
        xaxis_title='Position',
        yaxis_title='Reads',
        barmode='stack',  # Stack bars to show proportions
        dragmode='pan',  # Allows dragging to pan around
        xaxis_rangeslider_visible=True  # Adds a range slider at the bottom
    )
    
    # Show the figure
    fig.write_image(f'{output_dir}/Base_proportions.png')




# In[36]:


def regular_accuracy_plot(reference, result, output_dir):
    
    x_vals = []
    correct_score = []
    for ref_idx in range(0,len(reference)):
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
                    lw=0.5  # Line width
                )
            )
            theta-=math.pi/6
    plt.title('Position vs accuracy')
    plt.xlabel('Position')
    plt.ylabel('Proportion of bases at position that match refeernce sequence')
    plt.grid(True)
    plt.savefig(f'{output_dir}/Regular_accuracy.png')


# Create pie chart from tsv file
# 6 mutation types long, short: del, ins, wt, complex.
# 2 pie charts one mod, wt.one in depth

# In[6]:


def mutation_pie_chart(tsv_file, output_dir):

    # Initialize a list to store values from the second column
    second_column_values = []
    
    # Read the TSV file
    with open(tsv_file, mode='r') as file:
        tsv_reader = csv.reader(file, delimiter='\t')
        next(tsv_reader)  # Skip the header row
        for row in tsv_reader:
            # Assuming the second column is at index 1
            value = row[1]
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
    wt_count = value_counts.get('WT', 0)
    other_count = sum(count for key, count in value_counts.items() if key != 'WT')
    wt_vs_other = {'UNMODIFIED': wt_count, 'MODIFIED': other_count}
    
    # Dictionary for the second pie chart (each key in the dictionary)
    detailed_counts = dict(value_counts)
    
    # Create the first pie chart (WT vs Other)
    fig1 = go.Figure(data=[go.Pie(labels=list(wt_vs_other.keys()), values=list(wt_vs_other.values()))])
    fig1.update_layout(title='Proportion of WT vs Other')
    
    fig1.update_traces(
            textposition="outside",
            texttemplate="%{label} %{value} %{percent:.2p}",
    )
    
    
    # Create the second pie chart (detailed counts)
    fig2 = go.Figure(data=[go.Pie(labels=list(detailed_counts.keys()), values=list(detailed_counts.values()))])
    fig2.update_layout(title='Detailed Counts')
    
    fig2.update_traces(
            textposition="outside",
            texttemplate="%{label} %{value} %{percent:.2p}",
    )
    
    # Display the plots
    fig1.write_image(f'{output_dir}/Mutation_pie_chart.png')
    fig2.write_image(f'{output_dir}/Pattern_pie_chart.png')


# Small/Large Insertion pos graph

# In[7]:


def indel_per_position(tsv_file, reference, output_dir):

    
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['Read_id', 'Classification', 'Mutation_info'])
    
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
        total_reads+=1
        if row['Classification'] != 'WT':
            parsed = row['Mutation_info'].split(",")
            for mut in parsed:
                if "Del" in str(mut):
                    start_pos = int(re.match(r"(\d+)_", mut).group(1))
                    end_pos = int(re.search(r"_(\d+)", mut).group(1))
                    for position in range(start_pos, end_pos+1):
                        deletion_PerPos_dictionary[position] += 1
                if "Ins" in str(mut):
                    start_pos = int(re.match(r"(\d+)_", mut).group(1))
                    end_pos = int(re.search(r"_(\d+)", mut).group(1))
                    for position in range(start_pos, end_pos+1):
                        insertions_PerPos_dictionary[position] += 1
                if "Sub" in str(mut):
                    start_pos = int(re.match(r"(\d+)_", mut).group(1))
                    end_pos = int(re.search(r"_(\d+)", mut).group(1))
                    for position in range(start_pos, end_pos+1):
                        substitutions_PerPos_dictionary[position] += 1
    
    
    
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
                             line=dict(color='green'),  # Different color for insertions
                             name='Substitutions'))
    
    # Customize the layout
    fig.update_layout(
        title='Deletions and Insertions per Position',
        xaxis_title='Position',
        yaxis_title='Sequence %',
        showlegend=True,
        dragmode='pan',  # Allows dragging to pan around
        xaxis_rangeslider_visible=True  # Adds a range slider at the bottom
    )
    
    fig.update_yaxes(tickvals=y_ticks, ticktext=y_tick_labels)
    
    # Show the plot
    fig.write_image(f'{output_dir}/Deletion_and_insertions_per_position.png')


# Insertion average length vs position

# In[32]:


def Insertion_length(tsv_file, output_dir):
    
    # Read the data
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['Read_id', 'Classification', 'Mutation_info'])
    
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
                    insertions_PerPos_dictionary[start_pos] += 1
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
        margin=dict(l=0, r=0, t=50, b=50),  # Adjust margins if needed
        xaxis_rangeslider_visible=True  # Adds a range slider at the bottom
    )
    
    # Show the plot
    fig.write_image(f'{output_dir}/insertion_length.png')


# Del average length by pos

# In[33]:


def Deletion_length(tsv_file, output_dir):
    
    # Read the data
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['Read_id', 'Classification', 'Mutation_info'])
    
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
                    deletions_PerPos_dictionary[start_pos] += 1
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
        margin=dict(l=0, r=0, t=50, b=50),  # Adjust margins if needed
        xaxis_rangeslider_visible=True  # Adds a range slider at the bottom
    
    )
    
    # Show the plot
    fig.write_image(f'{output_dir}/deletion_length.png')


# Count of deletions of each size

# In[34]:


def Deletion_count_length(tsv_file, output_dir):

    # Read the data
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['Read_id', 'Classification', 'Mutation_info'])
    
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
                    length = end_pos - start_pos + 1
                    if length not in deletion_lengths_count:
                        deletion_lengths_count[length] = 0
                    deletion_lengths_count[length] += 1
    
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
        margin=dict(l=0, r=0, t=50, b=50),  # Adjust margins if needed
        xaxis_rangeslider_visible=True
    )
    
    # Show the plot
    fig.write_image(f'{output_dir}/Deletion_count_length.png')


# In[35]:


def insertion_count_length(tsv_file, output_dir):

    # Example data
    tsv_file = "/Users/benvyshedskiy/postproccesingpipeline/mutations_output_BCL11A_RNP2_T3.tsv"
    
    # Read the data
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['Read_id', 'Classification', 'Mutation_info'])
    
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
        xaxis_rangeslider_visible=True  # Adjust margins if needed
    )
    
    # Show the plot
    fig.write_image(f'{output_dir}/insertion_count_length.png')


# In[3]:


def allele_plot(ref_seq, cv_pos, input_file, output_dir):
        
    mut_dict = {}
    all_cnt = 0
    for line in open(input_file).readlines()[1:]:
        mut = line.strip().split()[2]
        all_cnt += 1
        if mut not in mut_dict:
            mut_dict[mut] = 1
        else:
            mut_dict[mut] += 1

    mut_list = []
    for i in sorted(mut_dict.items(), key=lambda x: x[1], reverse=True):
        mut_list.append([i[0], round(i[1]*100/all_cnt,2)])
    
    plot_window = 30
    show_line = 30
    
    ref_plot = ref_seq[cv_pos - plot_window: cv_pos + plot_window]
    
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
    
    fig, ax = plt.subplots(figsize=(12, 3 * show_line/10))
    
    
    #for info_n, info in enumerate(sorted(output_res.items(), key=lambda x: x[1], reverse=True)[:show_line]):
    for info_n, info in enumerate(mut_list[:show_line]):
        
        seq = ref_plot
        
        x_pos = 0.01
        y_pos = 0.9 - info_n * 0.1
        
        adjust_ins = 0
        print(info)
        
        if info[0] == 'None':
            for nucleotide in seq:
                ax.text(x_pos, y_pos, nucleotide, fontsize=14, fontfamily='monospace',
                        ha='left', va='center', bbox=dict(facecolor=nucleotide_colors[nucleotide], edgecolor='none', pad=2))
                x_pos += 0.0185 # Adjust spacing between nucleotides
            ax.text(0.0185 * plot_window * 2 +0.1, y_pos, f'{info[1]}% WT', fontsize=12, ha='left', va='center')
            continue
        
        mut_str = ''
        
        for mut in info[0].split(','):
    
            if mut.find('Sub') != -1:
                st_pos = int(mut.split(':')[0]) - (cv_pos - plot_window) + adjust_ins
                ed_pos = len(ref_seq)
            else:
                st_pos = int(mut.split('_')[0]) - (cv_pos - plot_window) + adjust_ins
                ed_pos = int(mut.split(':')[0].split('_')[1]) - (cv_pos - plot_window) + adjust_ins
    
            mut_type = mut.split(':')[1].replace('Large_','Large').split('_')[0]
            
            if st_pos > 60:
                if mut_type == 'Sub':
                    mut_str += f'{mut.split(":")[1].split("_")[1]}Sub,'
                else:
                    mut_len = int(mut.split(':')[1].replace('Large_','Large').split('_')[1])
                    if mut_type == 'Del':
                        mut_str += f'{mut_len}Del,'
                    elif mut_type == 'Ins':
                        mut_str += f'{mut_len}Ins,'
                    elif mut_type == 'LargeDel':
                        mut_str += f'{mut_len}Large_Del,'
                    elif mut_type == 'LargeIns':
                        mut_str += f'{mut_len}Large_Ins,'
                continue
            
            if st_pos < 0:
                st_pos = 0
            if ed_pos >= plot_window*2:
                ed_pos = plot_window*2 - 1
            
            if mut_type == 'Del':
                del_len = int(mut.split(':')[1].split('_')[1])
                seq = seq[:st_pos] + '-'*(ed_pos - st_pos + 1) + seq[ed_pos+1:]
                mut_str += f'{del_len}Del,'
            
            elif mut_type == 'Ins':
                ins_len = int(mut.split(':')[1].split('_')[1])
                rect = patches.Rectangle((0.008+st_pos*0.0185,  1-0.04-0.1*(info_n+1)), ins_len*0.0185, 0.1, linewidth=2, edgecolor='r', facecolor='none', zorder=10)
                ax.add_patch(rect)
                seq = (seq[:st_pos] + 'N' + seq[st_pos:])[:plot_window*2]
                mut_str += f'{ins_len}Ins,'
                adjust_ins += ins_len
            
            elif mut_type == 'LargeDel':
                del_len = int(mut.split(':')[1].replace('Large_','Large').split('_')[1])
                seq = seq[:st_pos] + ' '*(ed_pos - st_pos + 1) + seq[ed_pos+1:]
                ax.annotate('', xy=(0.005+st_pos*0.0185, 0.9-0.1*info_n), xytext=(0.01+(ed_pos+1)*0.0185, 0.9-0.1*info_n), arrowprops=dict(arrowstyle='-'), zorder=10)
                ax.text(0.005+(st_pos+ed_pos)/2*0.0185, 0.9-0.1*info_n, f'{del_len}bp', fontsize=10, backgroundcolor='white', ha='left', va='center', zorder=11,bbox=dict(facecolor='white', edgecolor='none', pad=2))
                mut_str += f'{del_len}Large_Del,'
                
            elif mut_type == 'LargeIns':
                ins_len = int(mut.split(':')[1].replace('Large_','Large').split('_')[1])
                seq = (seq[:st_pos] + 'N'*ins_len + seq[st_pos:])[:plot_window*2]
                if st_pos + ins_len > plot_window *2:
                    ins_len = plot_window*2 - st_pos
                rect = patches.Rectangle((0.008+st_pos*0.0185,  1-0.04-0.1*(info_n+1)), ins_len*0.0185, 0.1, linewidth=2, edgecolor='r', facecolor='none', zorder=10)
                ax.add_patch(rect)
                mut_str += f'{ins_len}Large_Ins,'
                adjust_ins += ins_len
            
            elif mut_type == 'Sub':
                mut_nt = mut.split('>')[1]
                seq = seq[:st_pos] + mut_nt + seq[st_pos+len(mut_nt):]
                print(seq)
                mut_str += f'{mut.split(":")[1].split("_")[1]}Sub,'
                
        for nucleotide in seq:
            ax.text(x_pos, y_pos, nucleotide, fontsize=14, fontfamily='monospace',
                    ha='left', va='center', bbox=dict(facecolor=nucleotide_colors[nucleotide], edgecolor='none', pad=2))
            x_pos += 0.0185 # Adjust spacing between nucleotides
        ax.text(0.0185 * plot_window * 2 +0.1, y_pos, f'{info[1]}% {mut_str[:-1]}', fontsize=12, ha='left', va='center')
    
    ax.annotate('', xy=(0.005+plot_window*0.0185, 1), xytext=(0.01+plot_window*0.0185, 1 - show_line/10 -0.1), arrowprops=dict(arrowstyle='-', linestyle='--', color='gray'), zorder=10)
    plt.xlim(0,0.005+(plot_window*2+1)*0.0185)
    plt.ylim(1 - show_line/10 -0.1, 1 )
    plt.axis('off')
    plt.savefig(output_dir + '/allel_plot.png')


# In[ ]:


def LD_tornado(tsv_file, cv_pos, ref_len, strand, output_dir):
        
    fig, ax = plt.subplots(figsize=(12, 10))
    
    LD_list =[]

    for line in open(tsv_file).readlines():
        line = line.split()
        if line[2].find('Large_Del') == -1:
            continue
        LD_st = int(line[2].split('_')[0])
        LD_ed = int(line[2].split(':')[0].split('_')[1])
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
    
    plt.savefig(output_dir + '/Large_deletion_tornado.png')

