def check_correctness():
    """
    This function checks the correctness of the provided answer by:
    1. Defining the genetic sequences and the standard genetic code.
    2. Translating each DNA sequence into an amino acid sequence.
    3. Identifying the type of mutation and its severity for each mutant.
    4. Determining the most severe mutation (the one most likely to create a non-functional protein).
    5. Comparing this scientific conclusion with the provided answer.
    """

    # 1. Define the data
    intact_seq = "ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGT"
    mutants = {
        "Mutant 1": "ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 2": "ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 3": "ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGT",
        "Mutant 4": "ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGT"
    }
    
    # Mapping from the question's options to the mutant names
    option_map = {
        "A": "Mutant 3",
        "B": "Mutant 1",
        "C": "Mutant 2",
        "D": "Mutant 4"
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer = "C"

    # Standard DNA to Amino Acid translation table
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'Stop', 'TAG':'Stop',
        'TGC':'C', 'TGT':'C', 'TGA':'Stop', 'TGG':'W',
    }

    # 2. Translation function
    def translate(seq):
        """Translates a DNA sequence into a protein sequence."""
        protein = ""
        codons = [seq[i:i+3] for i in range(0, len(seq) - (len(seq) % 3), 3)]
        for codon in codons:
            amino_acid = genetic_code.get(codon, 'X')
            if amino_acid == 'Stop':
                return protein, "Nonsense"
            protein += amino_acid
        return protein, "Full-length"

    # 3. Analyze each mutant
    results = {}
    intact_protein, _ = translate(intact_seq)

    for name, seq in mutants.items():
        protein, mutation_type = translate(seq)
        results[name] = {
            "protein": protein,
            "length": len(protein),
            "type": mutation_type,
            "description": ""
        }

    # Characterize each mutation
    # Mutant 1: CTC (Leu) -> TAC (Tyr) at codon 3. Missense.
    results["Mutant 1"]["description"] = f"Missense mutation. Protein length: {results['Mutant 1']['length']}."
    # Mutant 2: CTC (Leu) -> TAA (Stop) at codon 3. Nonsense.
    results["Mutant 2"]["description"] = f"Nonsense mutation (premature stop codon). Protein length: {results['Mutant 2']['length']}."
    # Mutant 3: Two missense mutations.
    results["Mutant 3"]["description"] = f"Multiple missense mutations. Protein length: {results['Mutant 3']['length']}."
    # Mutant 4: Missense + in-frame deletion.
    results["Mutant 4"]["description"] = f"Missense mutation and in-frame deletion. Protein length: {results['Mutant 4']['length']}."

    # 4. Determine the most severe mutation
    # A nonsense mutation that creates a severely truncated protein is the most severe.
    most_severe_mutant = None
    shortest_length = len(intact_protein)

    for name, data in results.items():
        if data["type"] == "Nonsense":
            if data["length"] < shortest_length:
                shortest_length = data["length"]
                most_severe_mutant = name
    
    if most_severe_mutant is None:
        return "Analysis Error: Could not identify a definitive knockout mutation among the choices."

    # 5. Compare with the provided answer
    correct_mutant_name = "Mutant 2"
    if most_severe_mutant != correct_mutant_name:
        return f"Code Logic Error: The analysis identified {most_severe_mutant} as the most severe, but it should be {correct_mutant_name}."

    llm_chosen_mutant = option_map.get(llm_answer)

    if llm_chosen_mutant is None:
        return f"Incorrect: The provided answer '{llm_answer}' is not a valid option (A, B, C, or D)."

    if llm_chosen_mutant == most_severe_mutant:
        return "Correct"
    else:
        return (f"Incorrect: The provided answer is '{llm_answer}', which corresponds to {llm_chosen_mutant}. "
                f"The correct choice is {most_severe_mutant}, which corresponds to option C. "
                f"Reason: {most_severe_mutant} is the only mutant with a nonsense mutation (a premature stop codon). "
                f"This mutation occurs at the 3rd codon, halting protein synthesis and producing a non-functional protein of only {shortest_length} amino acids. "
                f"This is the most reliable way to eliminate the protein's function.")

# Execute the check and print the result
print(check_correctness())