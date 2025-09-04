def check_answer():
    """
    This function checks the correctness of the given answer by simulating gene translation.
    The goal is to find the mutation that most probably inactivates a gene.
    The most effective way to do this is to introduce a premature stop codon (nonsense mutation).
    """

    # The provided answer from the LLM corresponds to Mutant 2
    llm_answer_choice = "D"
    answer_map = {"A": "Mutant 3", "B": "Mutant 1", "C": "Mutant 4", "D": "Mutant 2"}
    selected_mutant_name = answer_map.get(llm_answer_choice)

    # Define the DNA sequences for the intact gene and mutants
    sequences = {
        "Intact gene": "ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGT",
        "Mutant 1": "ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 2": "ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 3": "ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGT",
        "Mutant 4": "ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGT"
    }

    # Standard genetic code dictionary (DNA codon -> Amino Acid)
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'_STOP_', 'TAG':'_STOP_',
        'TGC':'C', 'TGT':'C', 'TGA':'_STOP_', 'TGG':'W',
    }

    # Function to translate a DNA sequence into a protein sequence
    def translate(dna_seq):
        protein = ""
        seq_len = len(dna_seq)
        for i in range(0, seq_len - (seq_len % 3), 3):
            codon = dna_seq[i:i+3]
            amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown codons
            if amino_acid == '_STOP_':
                protein += amino_acid
                break
            protein += amino_acid
        return protein

    # Find which mutant, if any, has a premature stop codon
    mutant_with_stop_codon = None
    intact_protein = translate(sequences["Intact gene"])
    
    for name, seq in sequences.items():
        if name == "Intact gene":
            continue
        
        protein = translate(seq)
        # A premature stop codon means the protein is shorter than the intact one and contains a stop signal
        if "_STOP_" in protein and len(protein) < len(intact_protein):
            mutant_with_stop_codon = name
            break

    # Check if the selected answer is the one with the premature stop codon
    if selected_mutant_name == mutant_with_stop_codon:
        return "Correct"
    else:
        reason = (f"The provided answer '{llm_answer_choice}' ({selected_mutant_name}) is incorrect.\n"
                  f"The question asks for the mutation most likely to eliminate the anti-nutritional compound, which requires inactivating the gene.\n"
                  f"Gene inactivation is most reliably achieved by a nonsense mutation that creates a premature stop codon, truncating the protein.\n"
                  f"Analysis of the sequences shows that Mutant 2's sequence ('ATGTTCTAAG...') translates to 'Met-Phe-STOP', introducing a stop codon at the third position.\n"
                  f"The selected answer, {selected_mutant_name}, does not cause this definitive inactivation. Therefore, the correct answer corresponds to Mutant 2, which is option D.")
        return reason

# Run the check and print the result
print(check_answer())