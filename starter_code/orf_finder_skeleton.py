# BIFS617 ORF Finder Skeleton
# Team members: Arjit, Seth, Tyler
# Date: 07/18/2024
# Fill in your code where prompted.

def load_fasta(filepath):
    # Team Member Name: Tyler
    # TODO: Parse a multi-line FASTA file, return dictionary {header: sequence}
    Pass
def load_fasta(filepath):
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


#Test if the code worked by printing the dictionary
fasta_dict = load_fasta("Example.fasta")
print(f"Loaded {len(fasta_dict)} sequences.")
print(fasta_dict)




def reverse_complement(seq):
    # Team Member Name: Arjit
    # TODO: Return reverse complement of sequence (optional: use Bio.Seq)
    

#First, I create a complement base map to get the reverse complement of whatever is loaded into the load_fasta function
Complement_map = {‘A’: ‘T’, ‘T’ : ‘A’, ‘C’: ‘G’, ‘G’:’C’, ‘N’:’N’}
# N is if a base is not ATC or G [usually cut off due to software or user error]
# Reverses the sequence using string slicing	
	reversed_seq = seq[::-1]
# Makes a list of complementary bases with all the adjustments
	complement_bases = [complement_map.get(base.upper(), ‘N’) for base in reversed_seq]
# Joins the list of bases back into a string (will be the complementary string)
	return “”.join(complement_bases)
Pass

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
                        else:           #Passes both: contains a stop codon and longer than minimum
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


def format_orf_output(header, frame, position, seq):
    # Team Member Name: Arjit
    # TODO: Return formatted FASTA header and codon-separated sequence

# This structures the important information in the header of the output
header = f”>{seq_id} | FRAME = {frame} | POS = {position} | LEN = {length} | {direction}”

# This adds the orf sequence in the proper place after the header
structured_ entry = f”{header}\n{seq}\n”

# This returns a properly formatted FASTA header and codon-seperated sequence
return structured_entry

pass

def create_visualization(orf_data, output_path):
    # Team Member Name: Tyler
    # TODO: create a visualization, save the file, for your ORF output
    pass

def main():
    # Team Member Name: Seth, Arjit, Tyler
    # TODO: Implement user input, sequence processing, and ORF printing and save the file
    stored_reverse_comp = {}
    stored_complete_header = []
    stored_complete_fra = {[]}
    stored_complete_pos = {[]}
    stored_complete_len = {[]}
    stored_complete_dir = {[]}
    stored_complete_ORF = {[]}
    
    user_minlen = input("Please enter a minimum length for ORFs:")
    user_filepath = input("Please enter a FASTA formatted file for ORF check:")
    
    with open(output_file, "w") as f:   #Clears file once when first running
        f.close()
        
    while user_filepath.lower() != "end" and user_filepath.lower() != "end":  #Runs until user inputs end characters
        try:    #Used to avoid breaking program if wrong file entered
            stored_fasta_info = load_fasta(user_filepath)  #headers and sequences stored from FASTA file in dictionary

            for key, value in stored_fasta_info:

                positive_complete_ORF = find_orfs(key, value, user_minlen)
                stored_complete_header = #Need to convert code to also return a list of header just like the rest and add here
                stored_complete_fra[key] = positive_complete_ORF[0]
                stored_complete_pos[key] = positive_complete_ORF[1]
                stored_complete_len[key] = positive_complete_ORF[2]
                stored_complete_dir[key] = positive_complete_ORF[3]
                stored_complete_ORF[key] = positive_complete_ORF[4]

                stored_reverse_comp[key] = reverse_complement(value)    #stores dictionary of repeat header, but with reverse complement
                negative_complete_ORF = find_orfs(key, stored_reverse_comp[key], user_minlen, "-")
                stored_complete_fra[key].append(negatvive_complete_ORF[0])
                stored_complete_pos[key].append(negatvive_complete_ORF[1])
                stored_complete_len[key].append(negatvive_complete_ORF[2])
                stored_complete_dir[key].append(negatvive_complete_ORF[3])
                stored_complete_ORF[key].append(negatvive_complete_ORF[4])

                user_filepath = input("ORF check complete. enter a new file or type 'End' to quit:")
                    
        except FileNotFoundError:
            user_file = input("File not found, enter a correct file name or type 'End' to quit:")
    pass


#Test if the code worked by printing the
fasta_dict = load_fasta("Example.fasta")
print(f"Loaded {len(fasta_dict)} sequences.")
print(fasta_dict)

if __name__ == "__main__":
    main()

