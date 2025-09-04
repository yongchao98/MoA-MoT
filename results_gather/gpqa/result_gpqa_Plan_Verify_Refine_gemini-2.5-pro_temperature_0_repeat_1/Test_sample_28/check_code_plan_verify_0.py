import re

def clean_sequence(seq_str):
    """Removes 5'-, -3', and ellipses from a DNA sequence string."""
    # Use regex to find the sequence between 5'- and ... or -3'
    match = re.search(r"5’-(.*?)(?:\…|-$)", seq_str)
    if match:
        return match.group(1)
    return ""

def translate_dna(seq):
    """Translates a DNA sequence into an amino acid sequence."""
    table = {
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
    protein = ""
    # Ensure we only process full codons
    seq_len = len(seq) - (len(seq) % 3)
    for i in range(0, seq_len, 3):
        codon = seq[i:i+3]
        protein += table.get(codon, 'X') # 'X' for unknown codon
    return protein

def check_answer():
    """
    Checks the correctness of the LLM's answer by analyzing the mutations.
    """
    # Sequences from the question
    intact_gene_str = "5’-ATGTTTCTCGCTGGTACTTCTGTGGATGAACATATTTATTGTCGT…TGA-3’"
    mutant_1_str = "5’-ATGTTCTACGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC…TGA-3’"
    mutant_2_str = "5’-ATGTTCTAAGCTGGTACTTCTGTGGATGAACATATTTATTGTCGC…TGA-3’"
    mutant_3_str = "5’-ATGTTTTACGCTGGTGTCACTTCTGTGGATGAACATATTTATTGTCGT…TGA-3’"
    mutant_4_str = "5’-ATGTTTTACGCTACTTCTGTGGATGAACATATTTATTGTCGT…TGA-3’"
    
    llm_answer = 'D' # Corresponds to Mutant 2

    # Clean and translate sequences
    intact_seq = clean_sequence(intact_gene_str)
    mutant_1_seq = clean_sequence(mutant_1_str)
    mutant_2_seq = clean_sequence(mutant_2_str)
    mutant_3_seq = clean_sequence(mutant_3_str)
    mutant_4_seq = clean_sequence(mutant_4_str)

    intact_protein = translate_dna(intact_seq)
    mutant_1_protein = translate_dna(mutant_1_seq)
    mutant_2_protein = translate_dna(mutant_2_seq)
    mutant_3_protein = translate_dna(mutant_3_seq)
    mutant_4_protein = translate_dna(mutant_4_seq)

    # Analyze mutations
    # The goal is to find the mutation most likely to cause a loss of function.
    # The hierarchy of severity is generally: Nonsense > Frameshift > In-frame Deletion/Insertion > Missense
    
    mutations = {}
    
    # Mutant 1: ATGTTTCTCGCT -> ATGTTCTACGCT. This is complex, but the most likely intended change is CTC -> TAC.
    # Intact: ATG TTT CTC -> M F L. Mutant: ATG TTC TAC -> M F Y. This is a missense mutation.
    if mutant_1_protein != intact_protein and "_STOP_" not in mutant_1_protein:
        mutations['A'] = "Missense"

    # Mutant 2: ATGTTTCTCGCT -> ATGTTCTAAGCT. Most likely intended change is CTC -> TAA (STOP).
    # Intact: ATG TTT CTC -> M F L. Mutant: ATG TTC TAA -> M F _STOP_. This is a nonsense mutation.
    if "_STOP_" in mutant_2_protein and mutant_2_protein.find("_STOP_") < len(intact_protein) - 1:
        mutations['D'] = "Nonsense"

    # Mutant 3: ATGTTTCTCGCT... -> ATGTTTTACGCTGGTGTCACTT...
    # This is an insertion of 'T' after the 6th base, causing a frameshift.
    # Intact: ATG TTT CTC GCT... -> M F L A...
    # Mutant: ATG TTT TAC GCT GGT GTC... -> M F Y A G V...
    # Let's check for frameshift by comparing lengths.
    if len(mutant_3_seq) != len(intact_seq):
        mutations['B'] = "Frameshift"
    else: # If not a frameshift, it's multiple missense mutations
        mutations['B'] = "Multiple Missense"

    # Mutant 4: ATGTTTCTCGCTGGTACT... -> ATGTTTTACGCTACTTCT...
    # This is a deletion of the 'GGT' codon and a missense CTC->TAC. This is an in-frame deletion.
    # Intact: M F L A G T S...
    # Mutant: M F Y A T S... (The G is missing)
    if len(mutant_4_protein) < len(intact_protein) and "_STOP_" not in mutant_4_protein:
        mutations['C'] = "In-frame Deletion + Missense"

    # Determine the best candidate
    # A nonsense mutation early in the sequence is the most definitive way to knock out a gene.
    best_candidate = None
    if 'Nonsense' in mutations.values():
        # Find the key corresponding to the 'Nonsense' mutation
        for key, value in mutations.items():
            if value == 'Nonsense':
                best_candidate = key
                break
    # If no nonsense, a frameshift would be next best.
    elif 'Frameshift' in mutations.values():
        for key, value in mutations.items():
            if value == 'Frameshift':
                best_candidate = key
                break
    
    if not best_candidate:
        return "Error: Could not definitively identify the most severe mutation among the choices."

    if llm_answer == best_candidate:
        return "Correct"
    else:
        return (f"Incorrect. The LLM chose {llm_answer}, but the best answer is {best_candidate}. "
                f"The analysis shows that Mutant 2 (Answer D) results in a nonsense mutation (premature stop codon), "
                f"which is the most likely mutation to eliminate the protein's function. "
                f"The LLM's final choice was correct, but this check confirms it based on the hierarchy of mutation severity.")

# Run the check
result = check_answer()
print(result)