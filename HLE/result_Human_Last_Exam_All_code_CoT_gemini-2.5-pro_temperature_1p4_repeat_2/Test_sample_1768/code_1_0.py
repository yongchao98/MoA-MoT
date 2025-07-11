def solve_orf():
    """
    Finds the correct reading frame and translates a DNA sequence from the middle of an ORF.
    """
    # The given nucleotide sequence
    dna_sequence = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard DNA codon table to one-letter amino acid codes
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

    # Iterate through the three possible forward reading frames (0, 1, 2)
    for frame in range(3):
        protein_sequence = []
        # Get the sequence for the current frame
        frame_sequence = dna_sequence[frame:]
        
        # Translate codons to amino acids in chunks of 3
        for i in range(0, len(frame_sequence), 3):
            codon = frame_sequence[i:i+3]
            
            # Ensure we have a full codon before translating
            if len(codon) == 3:
                amino_acid = genetic_code.get(codon, '?') # Use '?' for any unknown codons
                protein_sequence.append(amino_acid)

        # Join the list of amino acids into a string
        protein_string = "".join(protein_sequence)

        # An ORF in the middle of a gene should not contain stop codons ('*').
        # The first frame found without a stop codon is the correct one.
        if '*' not in protein_string:
            print(protein_string)
            # Stop searching once the correct frame is found
            return

solve_orf()