# BIFS617 ORF Finder Skeleton
# Team members: Arjit, Seth, Tyler
# Date: 07/18/2024
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
    lengths, frames = orf_data           # unpack the two lists
    idx = range(1, len(lengths)+1)
    #create figure
    plt.figure()
    bars = plt.bar(idx, lengths)
    
    # colour bars by frame
    colour_map = {1: "cornflowerblue", # Frame +1
                  2: "darkorange",     # Frame +2
                  3: "seagreen",       # Frame +3
                  4: "crimson",        # Frame -1
                  5: "orchid",         # Frame -2
                  6: "gold"}          # Frame -3
    
    for bar, frm in zip(bars, frames):
        bar.set_color(colour_map.get(frm, "grey"))
    #create labels and title
    plt.xlabel("ORF number")
    plt.ylabel("Length (nt)")
    plt.title("ORF lengths by reading frame")
    legend_handles = [Patch(color=c, label=f"Frame {f}")
                      for f, c in colour_map.items()]
    plt.legend(handles=legend_handles, title="Reading frame")
    #fix formating, save, and finish making graph
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

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
    stored_len_fra = [[], []]

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
                output = format_orf_output(h_key, frame, pos, length, direction, orf_seq)
                print(output)

                # Write to output file (append mode)
                with open(OUTPUT_FILE, "a") as f_out:
                    f_out.write(output + "\n")

            # Visualization data prep
            stored_len_fra[0].extend(stored_complete_len)
            stored_len_fra[1].extend(stored_complete_fra)

            OUTPUT_VIS_FILE = f"output/visualizations/{header}_orf_lengths.png"
            create_visualization(stored_len_fra, OUTPUT_VIS_FILE)

            user_filepath = input("ORF check complete. Enter a new file or type 'End' to quit: ")

        except FileNotFoundError:
            user_filepath = input("File not found. Enter a correct file name or type 'End' to quit: ")


if __name__ == "__main__":
    main()
