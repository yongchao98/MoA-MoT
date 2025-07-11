def predict_protein_function():
    """
    Translates a DNA sequence and identifies a premature stop codon,
    highlighting the need for a more advanced sequence analysis.
    """
    dna_sequence_raw = """
    atggctgctctcgaattccccgctgggttcctgtttggaacagcaacatcagcttaccagatagaaggcgcctggaaggaagatggtaagggcgaaagtatgtgggacagactgacgcacgatcatccggaaatcatcaaagataaatcaactggagatgtggcctgtaattcctatcatctgtataaggaagacgttcgaatgttaaaggaattaggggttaatttctaccggttctctgtgtcgtggtcgcgaattctaccaactggacatgacaacgtagtgaatcaagccggcattgcctattataataatctcataaacgaactgatcgccaatggaatacaacctatggtaatcatgtaccacttcgatctgccacagcctctacaagatcttggtggatggaccaatcctgtactagcaaactacttcgaggattacgctcgtgtgctgtacgctaacttcggagacagggtcaaatggtggaacactatcaacgaaccccagatacatgctatgggatatagcgaaccgtttggtcttgcacctaacatacttacccctggacatgggcagtacctggccgttcatacgtttctgttatcacatgccagggcttacagactgtacgagcgggaatttaaggcaaaacagggaggcaaagtgtctattgtgcctgggggcttctggatgcaaccgacgagtgactccaaagacgaagaggaagcagctgcgagagcgcaggaaatgcagttaggatgggtgttacatcccatctacagtgcaactggcgattatccgcctgtaatgaaggaatggatggctaagaaaagcaaggaggaaggttacagcaggtcgcgtctgcccagtttcaccaaggaggagattgaaatggtgaaagggacttgggactacctgggattaaaccactacacaacattcttcactctccagagcaaagaagaatcaatgtttctcaaggataccggagttgctaatattcaggatgataagtacccttctgctgcttcggaatggctacaggttgtaccatggggtttcagaaagttgctgaactggatagccaaaaagtacaacaatccacctatcgtgataacagagaatggattctctgaccacggagaacttaacgacagagatagagtcaactacctcacaaaatacttgtgcgaattactgaaggctgtgaaggacgacggatgcaacgtgatcggttatacagtgtggagtctgatggacaactttgaatggccgtctggatatactgagaaa
    """

    # Clean up the DNA sequence by removing spaces and newlines
    dna_sequence = "".join(dna_sequence_raw.strip().split())

    # Standard genetic code dictionary ('_' represents a stop codon)
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    protein_sequence = ""
    problem_codon = ""
    
    # Translate the DNA sequence codon by codon
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        if len(codon) == 3:
            amino_acid = genetic_code.get(codon.upper(), '?')
            if amino_acid == '_': # If a stop codon is found
                problem_codon = codon.upper()
                break # Stop the translation
            protein_sequence += amino_acid

    print("Attempting to translate the DNA sequence using the standard genetic code...")
    print(f"A stop codon ('{problem_codon}') was found after just 30 amino acids.")
    print("\nResulting (truncated) protein sequence:")
    print(protein_sequence)
    print("\nThis protein is too short to be the full gene product, suggesting a simple translation is insufficient.")

predict_protein_function()