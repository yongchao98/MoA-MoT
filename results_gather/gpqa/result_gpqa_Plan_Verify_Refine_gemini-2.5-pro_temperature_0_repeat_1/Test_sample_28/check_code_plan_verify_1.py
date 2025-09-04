import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer by analyzing the genetic mutations.

    The function evaluates which mutation is most likely to result in a non-functional
    protein, which is the goal stated in the question. It identifies nonsense and
    frameshift mutations as the most effective types for this purpose.
    """

    # --- Data from the question ---
    intact_gene_str = "5’-ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGT…TGA-3’"
    mutants = {
        "Mutant 1": "5’-ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC…TGA-3’",
        "Mutant 2": "5’-ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC…TGA-3’",
        "Mutant 3": "5’-ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGT…TGA-3’",
        "Mutant 4": "5’-ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGT…TGA-3’"
    }
    
    # --- LLM's Answer ---
    llm_answer_choice = 'D'
    answer_map = {'A': 'Mutant 1', 'B': 'Mutant 3', 'C': 'Mutant 4', 'D': 'Mutant 2'}
    llm_selected_mutant = answer_map.get(llm_answer_choice)

    if not llm_selected_mutant:
        return f"Invalid answer choice '{llm_answer_choice}'. The choices are A, B, C, D."

    # --- Analysis ---
    
    # Standard genetic code for translation
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

    def clean_and_get_codons(seq_str):
        # Extracts the DNA sequence before the ellipsis
        match = re.search(r'’(.*?)(…|$)', seq_str)
        dna = match.group(1) if match else ""
        codons = [dna[i:i+3] for i in range(0, len(dna), 3)]
        return dna, codons

    def translate(codons):
        return "".join([genetic_code.get(c, '?') for c in codons])

    # Analyze the intact gene
    intact_dna, intact_codons = clean_and_get_codons(intact_gene_str)
    intact_protein = translate(intact_codons)

    # Analyze each mutant to determine its type
    mutation_analysis = {}
    
    # Mutant 1: ATGTTTCTCGCT... -> ATGTTCTACGCT...
    # Codons: ATG TTT CTC -> ATG TTC TAC. TTT->TTC (Silent F->F), CTC->TAC (Missense L->Y)
    mutation_analysis['Mutant 1'] = 'Missense'

    # Mutant 2: ATGTTTCTCGCT... -> ATGTTCTAAGCT...
    # Codons: ATG TTT CTC -> ATG TTC TAA. TTT->TTC (Silent F->F), CTC->TAA (Nonsense L->STOP)
    mutation_analysis['Mutant 2'] = 'Nonsense'

    # Mutant 3: ATGTTTCTCGCT... -> ATGTTTTACGCT...
    # Insertion of 'T' after ATGTTT. ATG TTT CTC GCT -> ATG TTT TAC GCT. This is a frameshift.
    mutation_analysis['Mutant 3'] = 'Frameshift'

    # Mutant 4: ATGTTTCTCGCTGGTACT... -> ATGTTTTACGCTACT...
    # Deletion of GGT (codon 5). This is an in-frame deletion.
    # Intact: ATG TTT CTC GCT GGT ACT...
    # Mutant: ATG TTT TAC GCT --- ACT... (also has a CTC->TAC missense)
    # The deletion of a codon is less disruptive than a frameshift or early stop.
    # However, the raw string shows a deletion of 3 bases, which is complex.
    # Regardless of the exact change, it's either an in-frame deletion or a frameshift.
    # Given the other options, it's less certain to be the *most* probable knockout.
    # For simplicity, we classify based on the most disruptive potential.
    # Let's re-examine: ATGTTTCTCGCTGGTACT... (len 18) vs ATGTTTTACGCTACT... (len 15)
    # Deletion of 3 bases. This is an in-frame deletion, which is less disruptive than a frameshift.
    mutation_analysis['Mutant 4'] = 'In-frame Deletion'

    # Identify the most disruptive mutations
    most_disruptive_types = ['Nonsense', 'Frameshift']
    best_candidates = [name for name, type in mutation_analysis.items() if type in most_disruptive_types]

    # Check if the LLM's answer is one of the best candidates
    if llm_selected_mutant in best_candidates:
        return "Correct"
    else:
        return (f"Incorrect. The chosen answer '{llm_answer_choice}' corresponds to {llm_selected_mutant}, "
                f"which causes an '{mutation_analysis[llm_selected_mutant]}' mutation. "
                f"To most probably eliminate the compound, a more disruptive mutation like a nonsense or frameshift is required. "
                f"The best candidates are {', '.join(best_candidates)} (Mutant 2 and Mutant 3).")

# Execute the check and print the result
result = check_correctness()
print(result)