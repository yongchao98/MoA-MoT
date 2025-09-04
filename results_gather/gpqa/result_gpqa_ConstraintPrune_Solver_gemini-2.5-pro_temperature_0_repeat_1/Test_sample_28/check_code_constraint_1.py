import textwrap

def check_correctness():
    """
    Analyzes the mutations to determine which is most likely to knock out gene function
    and checks the correctness of the provided LLM answer and its reasoning.
    """
    # --- Data from the question ---
    intact_gene = "ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGT"
    mutants = {
        "Mutant 1": "ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 2": "ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 3": "ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGT",
        "Mutant 4": "ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGT"
    }
    llm_answer_choice = "C"
    option_map = {"A": "Mutant 1", "B": "Mutant 3", "C": "Mutant 2", "D": "Mutant 4"}

    # --- Genetic Code and Translation Function ---
    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_STOP_', 'TAG':'_STOP_', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GTC':'V', 'GTT':'V', 'GTA':'V', 'GTG':'V',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'AGA':'R', 'AGG':'R', 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'TGA':'_STOP_', 'TGC':'C', 'TGG':'W', 'TGT':'C', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G'
    }

    def translate(dna_seq):
        protein = []
        codons = textwrap.wrap(dna_seq, 3)
        for codon in codons:
            if len(codon) == 3:
                amino_acid = genetic_code.get(codon, 'X')
                protein.append(amino_acid)
                if amino_acid == '_STOP_':
                    break
        return "".join(protein)

    # --- Analysis Function ---
    def analyze_mutation(intact_seq, mutant_seq):
        len_diff = len(mutant_seq) - len(intact_seq)
        if len_diff != 0:
            if len_diff % 3 != 0:
                return "Frameshift", "Highly disruptive"
            else:
                return "In-frame indel", "Less disruptive"
        else:
            intact_protein = translate(intact_seq)
            mutant_protein = translate(mutant_seq)
            if "_STOP_" in mutant_protein and "_STOP_" not in intact_protein[:len(mutant_protein)]:
                return "Nonsense", "Highly disruptive"
            elif intact_protein != mutant_protein:
                return "Missense", "Less disruptive"
            else:
                return "Silent", "Not disruptive"

    # --- Perform Correct Analysis ---
    correct_analysis = {}
    for name, seq in mutants.items():
        correct_analysis[name] = analyze_mutation(intact_gene, seq)

    # --- Check LLM's Reasoning ---
    # The LLM's reasoning claims Mutant 3 is a frameshift. Let's verify.
    mutant3_analysis = correct_analysis["Mutant 3"]
    if mutant3_analysis[0] == "In-frame indel":
        llm_mutant3_claim = "Frameshift mutation"
        return (f"The provided answer's reasoning is incorrect. It classifies Mutant 3 as a '{llm_mutant3_claim}', "
                f"but the sequence has an insertion of 3 bases ({len(mutants['Mutant 3'])} vs {len(intact_gene)}), "
                f"which is an '{mutant3_analysis[0]}'. An in-frame insertion is considered '{mutant3_analysis[1]}' "
                f"and is not a frameshift. This is a critical error in the analysis, even though the final choice happens to be correct.")

    # --- Determine the Best Mutant ---
    best_mutant = None
    # Priority: Nonsense > Frameshift
    for name, (m_type, disruptive_level) in correct_analysis.items():
        if m_type == "Nonsense":
            best_mutant = name
            break
    if not best_mutant:
        for name, (m_type, disruptive_level) in correct_analysis.items():
            if m_type == "Frameshift":
                best_mutant = name
                break
    
    # If no nonsense or frameshift, the question is flawed, but we proceed.
    if not best_mutant:
        # Fallback for this specific question: find the first highly disruptive one
        for name, (m_type, disruptive_level) in correct_analysis.items():
            if disruptive_level == "Highly disruptive":
                best_mutant = name
                break

    # --- Final Check ---
    if option_map.get(llm_answer_choice) == best_mutant:
        return "Correct"
    else:
        return f"Incorrect. The best candidate is {best_mutant}, which contains a {correct_analysis[best_mutant][0].lower()} mutation. The provided answer chose {option_map.get(llm_answer_choice)}."

# Execute the check and print the result
result = check_correctness()
print(result)