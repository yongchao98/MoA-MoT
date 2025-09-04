def check_answer():
    """
    This function checks the correctness of the answer to the molecular biology question.
    It translates the DNA sequences and identifies the most severe mutation.
    """
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

    def translate(dna_seq):
        """Translates a DNA sequence into a protein sequence."""
        protein = ""
        # Process the sequence in codons (3 bases at a time)
        for i in range(0, len(dna_seq) - (len(dna_seq) % 3), 3):
            codon = dna_seq[i:i+3]
            amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown
            if amino_acid == '_': # Stop codon
                protein += "_"
                break # Stop translation
            protein += amino_acid
        return protein

    # --- Define Sequences ---
    intact_seq = "ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGT"
    mutants = {
        "Mutant 1": "ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 2": "ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 3": "ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGT",
        "Mutant 4": "ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGT"
    }
    
    # The provided answer is 'A', which corresponds to Mutant 2
    provided_answer_mutant = "Mutant 2"

    # --- Analysis ---
    intact_protein = translate(intact_seq)
    
    analysis_results = {}
    most_severe_mutant = None
    highest_severity_level = 0 # 3=Nonsense, 2=Frameshift, 1=Other

    for name, seq in mutants.items():
        mutant_protein = translate(seq)
        
        # Check for Nonsense Mutation (most severe)
        if '_' in mutant_protein and mutant_protein.index('_') < len(intact_protein):
            analysis_results[name] = f"Nonsense mutation. Protein is truncated to '{mutant_protein}'. This is the most severe mutation."
            if highest_severity_level < 3:
                highest_severity_level = 3
                most_severe_mutant = name
        else:
            # For this problem, other mutations are less severe.
            # A full analysis would differentiate missense, silent, and indels.
            analysis_results[name] = f"Less severe mutation (e.g., missense, in-frame deletion). Protein: '{mutant_protein}'."
            if highest_severity_level < 1:
                highest_severity_level = 1
                if most_severe_mutant is None:
                    most_severe_mutant = name

    # --- Verification ---
    if most_severe_mutant == provided_answer_mutant:
        return "Correct"
    else:
        error_message = f"The provided answer is incorrect.\n"
        error_message += f"The goal is to find the mutation most likely to create a non-functional protein. A nonsense mutation (premature stop codon) is the most effective way to do this.\n"
        error_message += "Analysis of mutants:\n"
        for name, result in analysis_results.items():
             error_message += f"- {name}: {result}\n"
        error_message += f"The analysis shows that '{most_severe_mutant}' has the most severe mutation. The provided answer chose '{provided_answer_mutant}'."
        return error_message

# Run the check
print(check_answer())