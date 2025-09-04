import textwrap

def check_correctness():
    """
    This function checks the correctness of the LLM's answer by analyzing the genetic mutations.
    It identifies the most disruptive mutation, which is the one most likely to result in a
    non-functional protein. The hierarchy of disruption is:
    1. Nonsense (premature stop codon) / Frameshift
    2. In-frame insertion/deletion
    3. Missense (amino acid change)
    4. Silent (no amino acid change)
    """

    # --- Data from the Question ---
    intact_gene = "ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGT"
    mutants = {
        "Mutant 1": "ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 2": "ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC",
        "Mutant 3": "ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGT",
        "Mutant 4": "ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGT"
    }
    llm_answer_choice = "C" # Corresponds to Mutant 2

    # --- Genetic Code and Translation Logic ---
    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'TCA':'S', 'TCC':'S',
        'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'CTA':'L', 'CTC':'L',
        'CTG':'L', 'CTT':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'_STOP_', 'TAG':'_STOP_', 'GAC':'D',
        'GAT':'D', 'GAA':'E', 'GAG':'E', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GTC':'V',
        'GTT':'V', 'GTA':'V', 'GTG':'V', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'AGA':'R',
        'AGG':'R', 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'TGA':'_STOP_', 'TGC':'C',
        'TGG':'W', 'TGT':'C', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G'
    }

    def translate(dna_seq):
        codons = textwrap.wrap(dna_seq, 3)
        protein = []
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
        
        # Frameshift (indel not multiple of 3)
        if len_diff != 0 and len_diff % 3 != 0:
            return "Frameshift mutation", 4
        
        # In-frame indel (indel is multiple of 3)
        if len_diff != 0 and len_diff % 3 == 0:
            return "In-frame insertion/deletion", 2

        # Substitutions
        intact_protein = translate(intact_seq)
        mutant_protein = translate(mutant_seq)
        
        # Nonsense (premature stop)
        if "_STOP_" in mutant_protein and mutant_protein.find("_STOP_") < len(intact_protein.replace("_STOP_", "")):
            return "Nonsense mutation", 5 # Highest disruption
        
        # Missense
        if intact_protein != mutant_protein:
            return "Missense mutation", 1
            
        # Silent
        return "Silent mutation", 0

    # --- Perform Analysis ---
    results = {}
    for name, seq in mutants.items():
        analysis, disruptiveness = analyze_mutation(intact_gene, seq)
        results[name] = {"analysis": analysis, "disruptiveness": disruptiveness}

    # Find the most disruptive mutation
    most_disruptive_mutant = max(results, key=lambda k: results[k]["disruptiveness"])

    # Map choice to mutant name
    choice_mapping = {"A": "Mutant 1", "B": "Mutant 3", "C": "Mutant 2", "D": "Mutant 4"}
    llm_answer_mutant = choice_mapping.get(llm_answer_choice)

    # --- Final Check ---
    if llm_answer_mutant == most_disruptive_mutant:
        return "Correct"
    else:
        reason = f"The provided answer '{llm_answer_choice}' ({llm_answer_mutant}) is incorrect.\n"
        reason += f"The most disruptive mutation is in {most_disruptive_mutant}, which has a '{results[most_disruptive_mutant]['analysis']}'.\n"
        reason += "This type of mutation is the most likely to result in a non-functional protein.\n"
        reason += "Analysis summary:\n"
        for name, res in results.items():
            reason += f"- {name}: {res['analysis']} (Disruption level: {res['disruptiveness']})\n"
        return reason

# Execute the check and print the result
print(check_correctness())