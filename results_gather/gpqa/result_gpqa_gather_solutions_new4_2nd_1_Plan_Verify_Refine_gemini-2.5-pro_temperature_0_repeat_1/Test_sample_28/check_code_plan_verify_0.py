import re

def check_answer():
    """
    Checks the correctness of the answer by analyzing the genetic mutations.
    The goal is to find the mutation most likely to create a non-functional protein.
    A nonsense mutation (premature stop codon) is the most effective way to achieve this.
    """
    
    # Standard DNA codon table
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
        # Remove non-alphabetic characters for clean processing
        seq = re.sub(r'[^ATGC]', '', seq.upper())
        protein = []
        # Ensure the sequence starts with a start codon
        if not seq.startswith('ATG'):
            return "No start codon"
        
        # Process in 3-base codons
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            if len(codon) == 3:
                amino_acid = genetic_code.get(codon, 'X') # 'X' for unknown codon
                protein.append(amino_acid)
                if amino_acid == '_STOP_':
                    break
        return "".join(protein)

    # Provided sequences
    intact_gene = "5’-ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGT…TGA-3’"
    mutant_1 = "5’-ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC…TGA-3’"
    mutant_2 = "5’-ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC…TGA-3’"
    mutant_3 = "5’-ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGT…TGA-3’"
    mutant_4 = "5’-ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGT…TGA-3’"
    
    sequences = {
        "Intact": intact_gene,
        "Mutant 1": mutant_1,
        "Mutant 2": mutant_2,
        "Mutant 3": mutant_3,
        "Mutant 4": mutant_4
    }

    # Translate all sequences
    translations = {name: translate_dna(seq) for name, seq in sequences.items()}
    
    # Identify the mutation type and severity
    analysis_results = {}
    best_mutant_name = None
    
    for i in range(1, 5):
        mutant_name = f"Mutant {i}"
        protein_seq = translations[mutant_name]
        
        # Check for nonsense mutation (premature stop codon)
        if "_STOP_" in protein_seq and len(protein_seq) < len(translations["Intact"].replace("_STOP_", "")):
            analysis_results[mutant_name] = f"Nonsense mutation. Protein: {protein_seq}. This introduces a premature stop codon, leading to a truncated, non-functional protein. This is the most effective knockout."
            best_mutant_name = mutant_name
            # A nonsense mutation is the most severe, so we can stop searching for the "best"
            break 
        else:
            # Simplified analysis for other types
            intact_dna_clean = re.sub(r'[^ATGC]', '', sequences["Intact"].upper())
            mutant_dna_clean = re.sub(r'[^ATGC]', '', sequences[mutant_name].upper())
            if len(mutant_dna_clean) < len(intact_dna_clean):
                 analysis_results[mutant_name] = f"In-frame deletion and/or missense mutation. Protein: {protein_seq}. Less likely to be a complete knockout than a nonsense mutation."
            else:
                 analysis_results[mutant_name] = f"Missense mutation(s). Protein: {protein_seq}. Less likely to be a complete knockout than a nonsense mutation."

    # The question maps options to mutants as follows:
    # A) Mutant 3, B) Mutant 4, C) Mutant 2, D) Mutant 1
    option_map = {'A': 'Mutant 3', 'B': 'Mutant 4', 'C': 'Mutant 2', 'D': 'Mutant 1'}
    
    # The provided answer is 'C'
    llm_answer_option = 'C'
    
    # Check if the LLM's choice corresponds to the best mutant found
    if best_mutant_name is None:
        return "Analysis Error: Could not definitively identify the most severe mutation. However, based on biological principles, a nonsense mutation should be the target."

    if option_map[llm_answer_option] == best_mutant_name:
        return "Correct"
    else:
        return (f"Incorrect. The analysis shows that {best_mutant_name} is the most likely to cause a gene knockout "
                f"because it results in a nonsense mutation (premature stop codon). "
                f"The LLM chose option {llm_answer_option}, which corresponds to {option_map[llm_answer_option]}, not {best_mutant_name}.")

# Run the check
result = check_answer()
print(result)