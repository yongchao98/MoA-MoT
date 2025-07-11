def translate_orf():
    """
    Finds the correct reading frame for a nucleotide sequence from an ORF
    and translates it into a protein sequence.
    """
    # The nucleotide sequence from the middle of an ORF
    dna_seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    
    # Standard genetic code mapping codons to one-letter amino acid codes
    gencode = {
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
    }

    # Iterate through the three possible forward reading frames (0, 1, 2)
    for frame_start_index in range(3):
        protein_seq = []
        has_stop_codon = False
        
        # Get the sequence for the current frame
        frame_dna = dna_seq[frame_start_index:].upper()
        
        # Process the sequence in chunks of 3 (codons)
        for i in range(0, len(frame_dna) - 2, 3):
            codon = frame_dna[i:i+3]
            amino_acid = gencode.get(codon, 'X') # 'X' for unknown codons
            
            # Check for a stop codon
            if amino_acid == '_':
                has_stop_codon = True
                break
            protein_seq.append(amino_acid)
        
        # If no stop codon was found, this is the correct frame
        if not has_stop_codon:
            final_protein = "".join(protein_seq)
            print(f"The correct reading frame is: {frame_start_index + 1}")
            print("The translated protein sequence is:")
            print(final_protein)
            # The puzzle has been solved, so we exit.
            return

    print("Could not find an open reading frame without a stop codon.")

# Run the function
translate_orf()
<<<HCCLAKTLNGTVKFSSSLLPYPNNLNSVLAERLEEGAGVG>>>