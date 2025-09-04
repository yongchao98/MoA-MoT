def check_correctness_of_llm_answer():
    """
    This function checks the correctness of the LLM's answer by analyzing the provided gene sequences.
    It identifies the most severe mutation and compares it with the LLM's choice.
    """

    # --- Data from the question ---
    intact_gene = "ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGT"
    mutants = {
        "Mutant 1": "ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 2": "ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 3": "ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGT",
        "Mutant 4": "ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGT",
    }

    # Mapping from the question's options to the mutant names
    options_map = {
        "A": "Mutant 3",
        "B": "Mutant 2",
        "C": "Mutant 1",
        "D": "Mutant 4",
    }

    # The final answer provided by the LLM to be checked
    llm_final_answer = "B"

    # --- Biological Analysis Logic ---
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

    def translate_dna(seq):
        """Translates a DNA sequence into a protein sequence."""
        protein = ""
        # Process in codons (3 bases) from the start
        for i in range(0, len(seq) - (len(seq) % 3), 3):
            codon = seq[i:i+3]
            amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown
            if amino_acid == '_STOP_':
                break
            protein += amino_acid
        return protein

    # The most severe mutation is a nonsense mutation (premature stop codon),
    # which results in a severely truncated protein. We can find this by
    # identifying the mutant that produces the shortest protein.
    
    intact_protein_len = len(translate_dna(intact_gene))
    shortest_protein_len = intact_protein_len
    most_severe_mutant_name = None
    
    for name, seq in mutants.items():
        protein = translate_dna(seq)
        if len(protein) < shortest_protein_len:
            shortest_protein_len = len(protein)
            most_severe_mutant_name = name

    # --- Verification ---
    # 1. Check if our analysis found a definitive knockout mutation.
    if most_severe_mutant_name is None:
        return "Analysis Error: No mutation was found to be more severe than others based on protein length. Expected a nonsense mutation."

    # 2. Check if the correct mutant was identified as the most severe.
    # Based on manual analysis, Mutant 2 has the nonsense mutation.
    if most_severe_mutant_name != "Mutant 2":
        return f"Analysis Error: The code identified {most_severe_mutant_name} as having the most severe mutation, but it should be Mutant 2."

    # 3. Check if the LLM's chosen option corresponds to the correct mutant.
    llm_chosen_mutant = options_map.get(llm_final_answer)
    if llm_chosen_mutant is None:
        return f"Invalid option '{llm_final_answer}' provided by the LLM. Options are A, B, C, D."
        
    if llm_chosen_mutant != most_severe_mutant_name:
        return f"Incorrect. The LLM chose option {llm_final_answer} which corresponds to {llm_chosen_mutant}. The correct choice is {most_severe_mutant_name} (Option B), which introduces a premature stop codon and is the most likely to eliminate the protein's function."

    # 4. Final check on the specific mutation in Mutant 2
    third_codon_mutant2 = mutants["Mutant 2"][6:9]
    if third_codon_mutant2 != "TAA":
        return f"Incorrect. The reasoning is flawed. The third codon in Mutant 2 should be 'TAA' (stop codon), but it is '{third_codon_mutant2}'."

    # If all checks pass, the LLM's answer and reasoning are correct.
    return "Correct"

# Run the check and print the result
print(check_correctness_of_llm_answer())