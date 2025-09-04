def check_answer():
    """
    Checks the correctness of the LLM's answer by simulating the biological impact of each mutation.
    The goal is to find the mutation most likely to create a non-functional protein.
    """

    # The LLM's answer choice to be verified
    llm_answer_choice = "C"

    # --- Data and Functions ---

    # Standard genetic code dictionary ('_' represents a stop codon)
    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'AGA':'R', 'AGG':'R',
        'GTC':'V', 'GTT':'V', 'GTA':'V', 'GTG':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    # Function to translate a DNA sequence into a protein sequence
    def translate(dna_seq):
        protein = ""
        seq_len = len(dna_seq)
        for i in range(0, seq_len - (seq_len % 3), 3):
            codon = dna_seq[i:i+3]
            amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown codons
            if amino_acid == '_': # Stop codon
                break
            protein += amino_acid
        return protein

    # DNA sequences from the question
    intact_gene = "ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGT"
    mutants = {
        1: "ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        2: "ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        3: "ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGT",
        4: "ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGT",
    }
    
    # Mapping from the answer choices to the mutant numbers
    choice_map = {"A": 1, "B": 3, "C": 2, "D": 4}

    # --- Analysis ---

    # Translate all sequences
    protein_intact = translate(intact_gene)
    proteins_mutant = {num: translate(seq) for num, seq in mutants.items()}

    # Identify the mutation types and their impact
    # Mutant 1: Missense mutation (changes an amino acid)
    # Mutant 2: Nonsense mutation (introduces a premature stop codon)
    # Mutant 3: In-frame insertion (adds an amino acid)
    # Mutant 4: In-frame deletion (removes an amino acid)

    # The most disruptive mutation is the one that causes premature termination (nonsense).
    # We can find this by identifying the mutant that produces the shortest protein.
    shortest_protein_len = len(protein_intact)
    most_disruptive_mutant_num = -1

    for num, protein in proteins_mutant.items():
        if len(protein) < shortest_protein_len:
            shortest_protein_len = len(protein)
            most_disruptive_mutant_num = num

    # --- Verification ---

    # Check 1: Was a most disruptive mutant found? It should be Mutant 2.
    if most_disruptive_mutant_num != 2:
        return f"Analysis Error: The code expected Mutant 2 to be the most disruptive (nonsense mutation), but instead found Mutant {most_disruptive_mutant_num}. Mutant 2 should result in a severely truncated protein."

    # Check 2: Does the LLM's answer correspond to the most disruptive mutant?
    selected_mutant_num = choice_map.get(llm_answer_choice)
    if selected_mutant_num is None:
        return f"Invalid Answer Format: The answer '{llm_answer_choice}' is not one of the valid choices (A, B, C, D)."

    if selected_mutant_num == most_disruptive_mutant_num:
        return "Correct"
    else:
        return (f"Incorrect. The answer identifies Mutant {selected_mutant_num} as the correct choice, "
                f"but the analysis shows that Mutant {most_disruptive_mutant_num} is the most likely to eliminate the compound. "
                f"This is because Mutant {most_disruptive_mutant_num} introduces a premature stop codon (a nonsense mutation), "
                f"which truncates the protein to only {shortest_protein_len} amino acids, almost certainly destroying its function. "
                f"The other mutations are less severe.")

# Run the check and print the result
print(check_answer())