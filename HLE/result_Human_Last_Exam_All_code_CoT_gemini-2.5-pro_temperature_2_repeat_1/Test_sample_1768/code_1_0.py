import textwrap

def find_and_translate_orf():
    """
    Finds the correct reading frame for a nucleotide sequence from the middle of an ORF
    and translates it into a protein sequence.
    """
    # The given nucleotide sequence
    dna_seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard DNA codon table
    genetic_code = {
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'TGC': 'C',
        'TGT': 'C', 'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'TTC': 'F', 'TTT': 'F', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G',
        'GGT': 'G', 'CAC': 'H', 'CAT': 'H', 'ATA': 'I', 'ATC': 'I',
        'ATT': 'I', 'AAA': 'K', 'AAG': 'K', 'CTA': 'L', 'CTC': 'L',
        'CTG': 'L', 'CTT': 'L', 'TTA': 'L', 'TTG': 'L', 'ATG': 'M',
        'AAC': 'N', 'AAT': 'N', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P',
        'CCT': 'P', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R',
        'CGG': 'R', 'CGT': 'R', 'AGA': 'R', 'AGG': 'R', 'TCA': 'S',
        'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'AGC': 'S', 'AGT': 'S',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'GTA': 'V',
        'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'TGG': 'W', 'TAC': 'Y',
        'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TGA': '*'
    }

    # Iterate through the three possible reading frames
    for frame in range(3):
        protein_sequence = ""
        codons = []
        # Slice the sequence for the current frame
        frame_seq = dna_seq[frame:]
        # Break the sequence into codons of 3 nucleotides
        for i in range(0, len(frame_seq) - 2, 3):
            codon = frame_seq[i:i+3]
            codons.append(codon)
            protein_sequence += genetic_code.get(codon.upper(), '?') # Use '?' for unknown codons

        # If a frame translates without stop codons, it's the correct one
        if '*' not in protein_sequence:
            print(f"Correct reading frame found: Frame {frame + 1}\n")
            
            print("Final Protein Sequence:")
            # Use textwrap to format the output neatly
            print(textwrap.fill(protein_sequence, width=60))
            print("\n" + "="*40)
            print("Translation Breakdown (Codon -> Amino Acid):")
            
            # Print the codon-by-codon translation
            translation_breakdown = [f"{codon} -> {genetic_code.get(codon.upper(), '?')}" for codon in codons]
            print("  ".join(translation_breakdown))
            
            # Once the correct frame is found and printed, exit the loop
            return

find_and_translate_orf()
<<<DIVVSEDLNGTVKFSSLLPYPNNLNSVLAERLEKW>>>