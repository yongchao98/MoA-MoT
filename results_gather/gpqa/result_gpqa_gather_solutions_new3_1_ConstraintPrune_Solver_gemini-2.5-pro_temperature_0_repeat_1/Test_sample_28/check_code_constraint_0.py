def check_answer():
    """
    Checks the correctness of the LLM's answer by analyzing the genetic mutations.
    """
    # Standard DNA codon table
    codon_table = {
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

    # Gene sequences from the question
    intact_gene = "ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGT"
    mutant_1 = "ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC"
    mutant_2 = "ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC"
    mutant_3 = "ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGT"
    mutant_4 = "ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGT"
    
    # The final answer provided by the LLM
    llm_answer = "C" # Corresponds to Mutant 2

    def translate(seq):
        """Translates a DNA sequence into an amino acid sequence."""
        protein = ""
        for i in range(0, len(seq) - (len(seq) % 3), 3):
            codon = seq[i:i+3]
            amino_acid = codon_table.get(codon, 'X') # 'X' for unknown codon
            if amino_acid == '_STOP_':
                protein += amino_acid
                break
            protein += amino_acid
        return protein

    # Analyze each sequence
    protein_intact = translate(intact_gene)
    protein_mutant_1 = translate(mutant_1)
    protein_mutant_2 = translate(mutant_2)
    protein_mutant_3 = translate(mutant_3)
    protein_mutant_4 = translate(mutant_4)

    analysis = {}
    
    # Mutant 1 Analysis: Missense mutation
    # Intact: ATG TTT CTC -> M F L
    # Mutant 1: ATG TTC TAC -> M F Y (CTC -> TAC is Leu -> Tyr)
    if 'Y' in protein_mutant_1[:3] and '_STOP_' not in protein_mutant_1:
        analysis['Mutant 1'] = "Missense mutation (Leucine -> Tyrosine). Less likely to be a complete knockout."
    
    # Mutant 2 Analysis: Nonsense mutation
    # Intact: ATG TTT CTC -> M F L
    # Mutant 2: ATG TTC TAA -> M F _STOP_ (CTC -> TAA is Leu -> STOP)
    if '_STOP_' in protein_mutant_2 and len(protein_mutant_2) < 5:
        analysis['Mutant 2'] = "Nonsense mutation (premature STOP codon). Most likely to cause complete loss of function."
        
    # Mutant 3 Analysis: Multiple missense mutations
    # Intact: ATG TTT CTC GCT GGT ACT -> M F L A G T
    # Mutant 3: ATG TTT TAC GCT GGT GTC -> M F Y A G V (CTC->TAC is Leu->Tyr; ACT->GTC is Thr->Val)
    if 'Y' in protein_mutant_3[:3] and 'V' in protein_mutant_3[:6] and '_STOP_' not in protein_mutant_3:
        analysis['Mutant 3'] = "Multiple missense mutations. Less certain to be a knockout than a nonsense mutation."

    # Mutant 4 Analysis: Missense + In-frame deletion
    # Intact: ATG TTT CTC GCT GGT ACT -> M F L A G T
    # Mutant 4: ATG TTT TAC GCT ACT -> M F Y A T (CTC->TAC is Leu->Tyr; GGT is deleted)
    # The resulting protein is shorter and has a substitution.
    if 'Y' in protein_mutant_4[:3] and protein_mutant_4[3:5] == protein_intact[3:4] + protein_intact[5:6]:
        analysis['Mutant 4'] = "Missense mutation plus an in-frame deletion. Less disruptive than a nonsense mutation."

    # Determine the most severe mutation
    most_severe_mutant = None
    for mutant, desc in analysis.items():
        if "Nonsense" in desc:
            most_severe_mutant = mutant
            break
    
    if not most_severe_mutant:
        return "Error in analysis: Could not definitively identify the most severe mutation."

    # Map mutant number to the option letter
    mutant_to_option = {"Mutant 1": "B", "Mutant 2": "C", "Mutant 3": "A", "Mutant 4": "D"}
    correct_option = mutant_to_option.get(most_severe_mutant)

    if llm_answer == correct_option:
        return "Correct"
    else:
        return (f"Incorrect. The LLM's answer is {llm_answer}, but the correct answer is {correct_option}.\n"
                f"Reason: The goal is to find the mutation most likely to create a non-functional protein. "
                f"The analysis shows that {most_severe_mutant} introduces a premature STOP codon (a nonsense mutation) "
                f"at the third codon position. This truncates the protein, making it non-functional. "
                f"This is the most severe and definitive mutation among the choices. The other mutations are missense "
                f"or in-frame deletions, which are less certain to cause a complete loss of function.")

# Run the check
print(check_answer())