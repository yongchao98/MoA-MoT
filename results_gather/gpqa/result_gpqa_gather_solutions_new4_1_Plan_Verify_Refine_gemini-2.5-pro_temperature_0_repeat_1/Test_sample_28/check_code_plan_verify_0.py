import re

def get_protein_sequence(dna_sequence):
    """Translates a DNA sequence into a protein sequence."""
    
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
    
    # Clean up the sequence by removing the ellipsis
    dna = re.sub(r'\.\.\..*', '', dna_sequence.upper())
    
    protein = ""
    # Find the start codon 'ATG'
    start_codon_index = dna.find('ATG')
    if start_codon_index == -1:
        return "No start codon found."
        
    # Start translation from the start codon
    for i in range(start_codon_index, len(dna), 3):
        codon = dna[i:i+3]
        if len(codon) == 3:
            amino_acid = codon_table.get(codon, 'X') # 'X' for unknown codon
            if amino_acid == '_STOP_':
                protein += amino_acid
                break
            protein += amino_acid
    return protein

def check_answer():
    """
    Checks the correctness of the provided answer by analyzing the mutations.
    """
    intact_gene = "5’-ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGT…TGA-3’"
    mutant_1 = "5’-ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC…TGA-3’"
    mutant_2 = "5’-ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC…TGA-3’"
    mutant_3 = "5’-ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGT…TGA-3’"
    mutant_4 = "5’-ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGT…TGA-3’"
    
    # The final answer provided by the LLM
    llm_answer = 'B' # Corresponds to Mutant 2

    # Translate all sequences
    p_intact = get_protein_sequence(intact_gene)
    p_mutant_1 = get_protein_sequence(mutant_1)
    p_mutant_2 = get_protein_sequence(mutant_2)
    p_mutant_3 = get_protein_sequence(mutant_3)
    p_mutant_4 = get_protein_sequence(mutant_4)

    analysis = []
    # Analyze Mutant 1
    if p_mutant_1 != p_intact and "_STOP_" not in p_mutant_1:
        analysis.append("Mutant 1: Missense mutation (changes amino acids but doesn't stop translation). Protein: " + p_mutant_1)
    
    # Analyze Mutant 2
    if "_STOP_" in p_mutant_2 and len(p_mutant_2) < len(p_intact):
        analysis.append("Mutant 2: Nonsense mutation (introduces a premature stop codon). Protein: " + p_mutant_2)
        most_severe_mutant = 2
    
    # Analyze Mutant 3
    if p_mutant_3 != p_intact and "_STOP_" not in p_mutant_3:
        analysis.append("Mutant 3: Missense mutations. Protein: " + p_mutant_3)

    # Analyze Mutant 4
    if p_mutant_4 != p_intact and "_STOP_" not in p_mutant_4:
        analysis.append("Mutant 4: Missense mutation and in-frame deletion. Protein: " + p_mutant_4)

    # Determine the correct choice based on biological principles
    # A nonsense mutation early in the sequence is the most likely to cause a complete loss of function.
    correct_choice_map = {1: 'D', 2: 'B', 3: 'C', 4: 'A'}
    correct_choice = correct_choice_map.get(most_severe_mutant)

    # Check if the LLM's answer matches the correct choice
    if llm_answer == correct_choice:
        return "Correct"
    else:
        reason = f"Incorrect. The provided answer is {llm_answer}, but the correct answer is {correct_choice}.\n"
        reason += "Reasoning:\n"
        for line in analysis:
            reason += f"- {line}\n"
        reason += "A nonsense mutation (Mutant 2) introduces a premature stop codon, leading to a truncated and non-functional protein. This is the most effective way to eliminate the gene's function among the given options."
        return reason

# Run the check
result = check_answer()
print(result)