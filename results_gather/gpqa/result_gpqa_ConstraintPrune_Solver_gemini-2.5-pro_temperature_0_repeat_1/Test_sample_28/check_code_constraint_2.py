import textwrap

def check_answer():
    """
    This function checks the correctness of the LLM's answer by analyzing the provided gene mutations.
    It determines which mutation is most likely to result in a non-functional protein.
    """

    # --- Data and Helper Functions ---

    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_STOP_', 'TAG':'_STOP_',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'AGA':'R', 'AGG':'R',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TGA':'_STOP_', 'TGC':'C', 'TGG':'W', 'TGT':'C',
        'GTC':'V', 'GTT':'V', 'GTA':'V', 'GTG':'V'
    }

    def translate(seq):
        """Translates a DNA sequence into an amino acid sequence, stopping at the first STOP codon."""
        seq = seq[:len(seq) - (len(seq) % 3)] # Ensure sequence is multiple of 3
        codons = textwrap.wrap(seq, 3)
        protein = []
        for codon in codons:
            amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown codons
            if amino_acid == '_STOP_':
                protein.append(amino_acid)
                break
            protein.append(amino_acid)
        return "".join(protein)

    def analyze_mutation(intact_seq, mutant_seq):
        """
        Analyzes the type and impact of a mutation.
        Returns a tuple: (mutation_type_string, impact_level_string)
        """
        len_diff = len(mutant_seq) - len(intact_seq)

        # 1. Check for Frameshift (indel not divisible by 3)
        if len_diff != 0 and len_diff % 3 != 0:
            return "Frameshift", "Highly disruptive"

        # 2. Check for In-frame indel (indel divisible by 3)
        if len_diff % 3 == 0 and len_diff != 0:
            m_type = "In-frame insertion" if len_diff > 0 else "In-frame deletion"
            return m_type, "Less disruptive"

        # 3. If length is the same, check for substitution effects
        protein_intact = translate(intact_seq)
        protein_mutant = translate(mutant_seq)

        # Check for Nonsense mutation (premature stop)
        if '_STOP_' in protein_mutant and len(protein_mutant) < len(protein_intact):
            return "Nonsense", "Highly disruptive"

        # Check for Missense/Silent mutation
        if protein_intact != protein_mutant:
            return "Missense", "Less disruptive"
        
        return "Silent", "Less disruptive"

    # --- Main Checking Logic ---

    # Sequences from the question
    intact_gene = "ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGT"
    mutants = {
        "Mutant 1": "ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 2": "ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 3": "ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGT",
        "Mutant 4": "ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGT"
    }
    
    # The LLM's final answer choice
    llm_final_choice = "C"
    option_map = {"A": "Mutant 1", "B": "Mutant 3", "C": "Mutant 2", "D": "Mutant 4"}
    llm_chosen_mutant = option_map.get(llm_final_choice)

    # Perform code-based analysis
    analysis_results = {}
    highly_disruptive_mutants = []
    
    for name, seq in mutants.items():
        m_type, impact = analyze_mutation(intact_gene, seq)
        analysis_results[name] = {"type": m_type, "impact": impact}
        if impact == "Highly disruptive":
            highly_disruptive_mutants.append(name)

    # --- Verification ---
    
    # Expected analysis based on biological principles:
    # Mutant 1: Missense (Less disruptive)
    # Mutant 2: Nonsense (Highly disruptive)
    # Mutant 3: In-frame insertion of 3bp (Less disruptive)
    # Mutant 4: In-frame deletion of 3bp (Less disruptive)
    
    # Check if the code's analysis matches expectations
    if (analysis_results["Mutant 1"]["impact"] != "Less disruptive" or
        analysis_results["Mutant 2"]["impact"] != "Highly disruptive" or
        analysis_results["Mutant 3"]["impact"] != "Less disruptive" or
        analysis_results["Mutant 4"]["impact"] != "Less disruptive"):
        return f"Reasoning Error: The code's analysis of mutation impacts does not match biological expectations. Results: {analysis_results}"

    # Check if there is exactly one highly disruptive mutant as expected
    if len(highly_disruptive_mutants) != 1:
        return f"Reasoning Error: Expected one highly disruptive mutation, but found {len(highly_disruptive_mutants)}: {highly_disruptive_mutants}"

    # Check if the LLM's choice matches the most disruptive mutant
    most_disruptive_mutant = highly_disruptive_mutants[0]
    if most_disruptive_mutant != llm_chosen_mutant:
        return (f"Incorrect final answer: The most disruptive mutation is {most_disruptive_mutant}, "
                f"but the LLM chose {llm_chosen_mutant} (Option {llm_final_choice}).")

    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)