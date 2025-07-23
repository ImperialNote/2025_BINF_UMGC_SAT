# BIFS617 ORF Finder Skeleton
# Team members: Arjit, Seth Nelson
# Date: 
# Fill in your code where prompted.

def load_fasta(filepath):
    # Team Member Name:
    # TODO: Parse a multi-line FASTA file, return dictionary {header: sequence}
    pass

def reverse_complement(seq):
    # Team Member Name:
    # TODO: Return reverse complement of sequence (optional: use Bio.Seq)
    pass

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
                            break
                    else:
                        continue
    
        z += 1
        r_frame += 1

    return frame_diction, pos_diction, len_diction, ORF_diction


def format_orf_output(header, frame, position, seq):
    # Team Member Name:
    # TODO: Return formatted FASTA header and codon-separated sequence
    pass

def create_visualization(orf_data, output_path):
    # Team Member Name:
    # TODO: create a visualization, save the file, for your ORF output
    pass

def main():
    # Team Member Name:
    # TODO: Implement user input, sequence processing, and ORF printing and save the file
    pass

if __name__ == "__main__":
    main()
