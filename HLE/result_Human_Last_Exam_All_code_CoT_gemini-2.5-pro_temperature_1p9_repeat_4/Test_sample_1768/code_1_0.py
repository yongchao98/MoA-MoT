import sys

def find_and_translate_orf():
    """
    Finds the correct reading frame and translates a nucleotide sequence into a protein sequence.

    The script iterates through all three reading frames, translates them, and identifies
    the correct frame as the one that does not contain a stop codon, as the sequence is
    from the middle of an ORF.
    """
    dna_seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    dna_upper = dna_seq.upper()

    # Standard genetic code (DNA codons to one-letter amino acid code)
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

    # Iterate through the three possible reading frames
    for frame in range(3):
        protein_seq = []
        has_stop_codon = False
        
        # Translate the sequence for the current frame
        for i in range(frame, len(dna_upper) - 2, 3):
            codon = dna_upper[i:i+3]
            amino_acid = codon_table.get(codon, '?') # '?' for unknown codons
            
            if amino_acid == '*':
                has_stop_codon = True
                break
            protein_seq.append(amino_acid)

        # An ORF from the middle of a gene should not have a stop codon.
        if not has_stop_codon:
            final_protein = "".join(protein_seq)
            print(f"Correct reading frame: {frame + 1}")
            print(f"Translated protein sequence:")
            print(final_protein)
            # Once the correct frame is found, we can stop.
            return

    # This part would be reached if no open reading frame is found
    print("No open reading frame without a stop codon was found.", file=sys.stderr)


# Execute the function
find_and_translate_orf()
<<<IVVSEDLNGTVKFSSSLLPYPNNLNSSVLAERLEKW>>>