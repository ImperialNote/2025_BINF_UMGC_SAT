# BIFS617 ORF Finder Skeleton
# Team members: Arjit, Seth, Tyler
# Date: 08/06/2024
# Fill in your code where prompted.

import matplotlib.pyplot as plt
import os
from matplotlib.patches import Patch

os.makedirs("output/orfs", exist_ok=True)
os.makedirs("output/visualizations", exist_ok=True)

OUTPUT_FILE = "output/orfs/orf_output.fasta"

with open(OUTPUT_FILE, "w") as f:   #Clears file once when first running
    f.close()

def load_fasta(filepath):
    # Team Member Name: Tyler
    # TODO: Parse a multi-line FASTA file, return dictionary {header: sequence}
    # We will store the results in the dictionary called sequences: {header: sequence}
    sequences = {}
    # Open the file for reading
    file = open(filepath, 'r')
    # Track the current header (our Key) and its sequence (the value), these start empty but will be filled up as we parse the fasta
    header = None
    sequence = ""
    # Read the file line by line
    for line in file:
        line = line.strip()           # remove spaces
        #  Header line
        if line.startswith('>'):
            # If we were collecting a previous record, save it
            if header is not None:
                sequences[header] = sequence
            # drop the leading '>'
            header = line[1:]
            # reset sequence holder        
            sequence = ""            
        # Sequence line
        elif line != "":              
            sequence = sequence + line
    # Save the final record after the loop ends
    if header != None:
        sequences[header] = sequence
    # Close the file and hand back the dictionary
    file.close()
    return sequences

def reverse_complement(seq):
    # Team Member Name: Arjit
    # TODO: Return reverse complement of sequence (optional: use Bio.Seq)
    
    #First, I create a complement base map to get the reverse complement of whatever is loaded into the load_fasta function
    Complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    # N is if a base is not ATC or G [usually cut off due to software or user error]
    # Reverses the sequence using string slicing	
    reversed_seq = seq[::-1]
    # Makes a list of complementary bases with all the adjustments
    complement_bases = [Complement_map.get(base.upper(), 'N') for base in reversed_seq]
    # Joins the list of bases back into a string (will be the complementary string)
    return ''.join(complement_bases)

def find_orfs(header, sequence, min_len, strand="+"):
    # Team Member Name: Seth    
    # TODO: Identify ORFs in all 3 reading frames for one strand

    #Establishes two arrays of codons and dictionaries for each value stored, in order and tied to header
    start_codons = ["ATG", "TTG", "GTG"]
    stop_codons = ["TAA", "TAG", "TGA"]
    ORF_diction = {header: []}
    frame_diction = {header: []}
    pos_diction = {header: []}
    len_diction = {header: []}
    dir_diction = {header: []}

    #Values used to iterrate through both text and while loop to truncate to one loop for all 3 reading frames
    z = 0
    if strand == "-":   #Checks if the the reading frames need to iterate through 4, 5, and 6 on negative strand
        r_frame = 4
    else:
        r_frame = 1
    #While loop which controls and iterates through the reading frame
    while z < 3:
        for i in range(z, len(sequence), 3):
            total_ORF = ""  #Reusable variable
            ORF_1_codon = sequence[i:i+3]
            if ORF_1_codon in start_codons: #If a codon every 3 bases matches start, will execute loop
                for p in range(i, len(sequence), 3):    #New range this time starting at the start codon
                    current_codon = sequence[p:p+3]     #Needs to stay in groups of 3 to continue checking for stop codons
                    total_ORF += current_codon          #Adds current 3 bases to running total after start codon discovery
                    if current_codon in stop_codons:    #Begins check to see if ORF meets minimum requirements
                        if len(total_ORF) <= min_len:        #Breaks loop if too short
                            break
                        else:           # es both: contains a stop codon and is longer than the minimum
                            ORF_diction[header].append(total_ORF)                           #Stores ORF and all other values in appropriate dictionaries
                            frame_diction[header].append(r_frame)
                            pos_diction[header].append(i + 1)
                            len_diction[header].append(len(total_ORF))
                            dir_diction[header].append(strand)
                            break
                    else:
                        continue
        z += 1
        r_frame += 1
    return [frame_diction, pos_diction, len_diction, dir_diction, ORF_diction]

def format_orf_output(header, frame, position, length, direction, seq):
    # Team Member Name: Arjit
    # TODO: Return formatted FASTA header and codon-separated sequence
    # Convert direction symbol to text
    if direction == "+":
        direction_str = "FOR"
    elif direction == "-":
        direction_str = "REV"
    else:
        direction_str = direction  # Fall back to whatever was passed if not +/- 
    # This structures the important information in the header of the output
    fasta_header = f">{header} | {frame} | {position} | {length} | {direction_str}"
    codon_seq = ' '.join(seq[i:i+3] for i in range(0, len(seq), 3))
    # This adds the orf sequence in the proper place after the header
    structured_entry = f"{fasta_header}\n{codon_seq}\n"
    # This returns a properly formatted FASTA header and codon-seperated sequence
    return structured_entry


def create_visualization(orf_data, output_path):
    #This will unpack all three lists of data
    all_lengths, all_frames, all_headers = orf_data
    if not all_lengths:
        print("No ORF data to visualize.")
        return
    # Group the flat data lists into a dictionary for easier handling
    grouped_data = {}
    unique_headers_in_order = []
    for header, length, frame in zip(all_headers, all_lengths, all_frames):
        if header not in grouped_data:
            grouped_data[header] = {'lengths': [], 'frames': []}
            unique_headers_in_order.append(header) # Preserve original order
        grouped_data[header]['lengths'].append(length)
        grouped_data[header]['frames'].append(frame)
    # Calculate bar positions, colors, and labels for the grouped plot
    plt.figure(figsize=(16, 8))
    colour_map = {1: "cornflowerblue",
                  2: "darkorange",
                  3: "seagreen",
                  4: "crimson",
                  5: "orchid",
                  6: "gold"}
    legend_handles = [Patch(color=c, label=f"Frame {f}") for f, c in colour_map.items()]
    all_x_positions = []
    all_plot_lengths = []
    all_colors = []
    group_tick_positions = []
    group_tick_labels = []
    current_pos = 1
    GROUP_GAP = 2 # The space between each sequence's bars
    for header in unique_headers_in_order:
        data = grouped_data[header]
        group_start_pos = current_pos     
        for i in range(len(data['lengths'])):
            all_x_positions.append(current_pos)
            all_plot_lengths.append(data['lengths'][i])
            all_colors.append(colour_map.get(data['frames'][i], 'grey'))
            current_pos += 1  
        group_end_pos = current_pos - 1
        tick_pos = (group_start_pos + group_end_pos) / 2.0
        group_tick_positions.append(tick_pos)
        group_tick_labels.append(header)        
        current_pos += GROUP_GAP
    #Plot the data and format the chart
    plt.bar(all_x_positions, all_plot_lengths, color=all_colors)
    plt.xlabel("Sequence Header")
    plt.ylabel("Length (nt)")
    plt.title("ORF Lengths by Sequence and Reading Frame")  
    ax = plt.gca()
    ax.set_xticks(group_tick_positions)
    ax.set_xticklabels(group_tick_labels, rotation=45, ha="right")
    ax.xaxis.set_minor_locator(plt.NullLocator())
    plt.legend(handles=legend_handles, title="Reading frame")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"\nGrouped visualization saved to {output_path}")


def main():
    # Team Member Name: Seth, Arjit, Tyler
    # TODO: Implement user input, sequence processing, and ORF printing and save the file

    stored_reverse_comp = {}
    stored_complete_fra = []
    stored_complete_pos = []
    stored_complete_len = []
    stored_complete_dir = []
    stored_complete_ORF = []
    stored_complete_headers = []

    # this is for visualization function. Uses input of two lists: lengths and frames
    vis_orf_data = [[], [], []]

    # default user_minlen set to 5 as instructed in project outline
    user_minlen = 5
    user_minlen = int(input("Please enter a minimum length for ORFs: "))
    user_filepath = input("Please enter a FASTA formatted file for ORF check: ")

    while user_filepath.lower() != "end":  # Runs until user inputs end
        try:  # Used to avoid breaking program if wrong file entered
            stored_fasta_info = load_fasta(user_filepath)  # headers and sequences stored from FASTA file

            for h_key, value in stored_fasta_info.items():

                #Positive strand ORFs
                positive_complete_ORF = find_orfs(h_key, value, user_minlen)
                stored_complete_fra.extend(positive_complete_ORF[0][h_key])
                stored_complete_pos.extend(positive_complete_ORF[1][h_key])
                stored_complete_len.extend(positive_complete_ORF[2][h_key])
                stored_complete_dir.extend(positive_complete_ORF[3][h_key])
                stored_complete_ORF.extend(positive_complete_ORF[4][h_key])
                stored_complete_headers.extend([h_key] * len(positive_complete_ORF[4][h_key]))

                #Negative strand ORFs
                stored_reverse_comp[h_key] = reverse_complement(value)
                negative_complete_ORF = find_orfs(h_key, stored_reverse_comp[h_key], user_minlen, "-")
                stored_complete_fra.extend(negative_complete_ORF[0][h_key])
                stored_complete_pos.extend(negative_complete_ORF[1][h_key])
                stored_complete_len.extend(negative_complete_ORF[2][h_key])
                stored_complete_dir.extend(negative_complete_ORF[3][h_key])
                stored_complete_ORF.extend(negative_complete_ORF[4][h_key])
                stored_complete_headers.extend([h_key] * len(negative_complete_ORF[4][h_key]))

            #Iterate and print/save all ORFs
            for i in range(len(stored_complete_ORF)):
                orf_seq = stored_complete_ORF[i]
                frame = stored_complete_fra[i]
                pos = stored_complete_pos[i]
                length = stored_complete_len[i]
                direction = stored_complete_dir[i]
                header = stored_complete_headers[i]
                
                #Stores output variable for both printing to screen and to file
                output = format_orf_output(header, frame, pos, length, direction, orf_seq)
                print(output)

                # Write to output file (append mode)
                with open(OUTPUT_FILE, "a") as f_out:
                    f_out.write(output + "\n")

            # Visualization data prep
            vis_orf_data[0].extend(stored_complete_len)
            vis_orf_data[1].extend(stored_complete_fra)
            vis_orf_data[2].extend(stored_complete_headers)

            OUTPUT_VIS_FILE = f"output/visualizations/{header}_orf_lengths.png"
            create_visualization(vis_orf_data, OUTPUT_VIS_FILE)

            user_filepath = input("ORF check complete. Enter a new file or type 'End' to quit: ")

        except FileNotFoundError:
            user_filepath = input("File not found. Enter a correct file name or type 'End' to quit: ")


if __name__ == "__main__":
    main()
