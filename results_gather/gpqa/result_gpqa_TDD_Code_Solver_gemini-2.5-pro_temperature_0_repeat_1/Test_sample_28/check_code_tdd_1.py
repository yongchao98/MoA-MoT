def check_mutation_answer():
    """
    Checks the correctness of the LLM's answer by analyzing the impact of each mutation.
    The function identifies the mutation most likely to create a non-functional protein,
    which is typically a nonsense mutation (premature stop codon).
    """
    # Standard DNA to Amino Acid translation table
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

    def translate(seq):
        """Translates a DNA sequence into an amino acid sequence."""
        # Ensure the sequence is a multiple of 3 for codon processing
        seq_len = len(seq)
        protein = []
        for i in range(0, seq_len - (seq_len % 3), 3):
            codon = seq[i:i+3]
            amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown codon
            protein.append(amino_acid)
            if amino_acid == '*':
                break  # Stop translation at a stop codon
        return "".join(protein)

    # --- Data from the question ---
    intact_gene = "ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGTTGA"
    mutants = {
        "Mutant 1": "ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGCTGA",
        "Mutant 2": "ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGCTGA",
        "Mutant 3": "ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGTTGA",
        "Mutant 4": "ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGTTGA",
    }
    llm_answer_option = "B"
    llm_answer_mutant = "Mutant 2"

    # --- Analysis ---
    # 1. Translate the intact gene to establish a baseline
    intact_protein = translate(intact_gene)
    # The length of the coding sequence before the natural stop codon
    intact_coding_length = intact_protein.find('*')
    if intact_coding_length == -1:
        intact_coding_length = len(intact_protein)

    # 2. Analyze each mutant to find the most severe mutation
    premature_stop_mutants = {}
    for name, seq in mutants.items():
        mutant_protein = translate(seq)
        stop_pos = mutant_protein.find('*')
        
        # Check if a premature stop codon exists
        if stop_pos != -1 and stop_pos < intact_coding_length:
            premature_stop_mutants[name] = stop_pos

    # 3. Determine the best candidate
    # A premature stop codon is the most likely knockout mutation.
    # If multiple exist, the earliest one is the most severe.
    if not premature_stop_mutants:
        return "Incorrect. The provided answer claims a nonsense mutation, but the code did not find any mutant with a premature stop codon."

    # Find the mutant with the earliest stop codon
    most_severe_mutant = min(premature_stop_mutants, key=premature_stop_mutants.get)

    # 4. Verify the LLM's answer
    if most_severe_mutant == llm_answer_mutant:
        # Let's double-check the specific codon change mentioned by the LLM
        # Intact gene 3rd codon: CTC (Leu)
        # Mutant 2 3rd codon: TAA (Stop)
        intact_codon_3 = intact_gene[6:9]
        mutant2_codon_3 = mutants["Mutant 2"][6:9]
        
        if intact_codon_3 == "CTC" and mutant2_codon_3 == "TAA":
            return "Correct"
        else:
            return f"Incorrect. The reasoning is flawed. The code identified {most_severe_mutant} as correct, but the specific codon change mentioned in the LLM's explanation is wrong. Expected CTC->TAA, but found {intact_codon_3}->{mutant2_codon_3}."
            
    else:
        return f"Incorrect. The analysis shows that '{most_severe_mutant}' is the most likely candidate to eliminate the compound because it has the earliest premature stop codon. The provided answer was '{llm_answer_mutant}' (Option {llm_answer_option})."

# Run the check
result = check_mutation_answer()
print(result)