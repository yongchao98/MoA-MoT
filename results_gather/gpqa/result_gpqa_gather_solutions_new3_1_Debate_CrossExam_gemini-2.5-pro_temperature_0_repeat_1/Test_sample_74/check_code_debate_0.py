import textwrap

def check_answer():
    """
    Checks the correctness of the LLM's answer by translating the provided DNA sequence.
    """
    # The DNA sequence from the problem description
    dna_sequence = "ATGTACCCATACGATGTTCCAGATTACGCCAAATGACTCTGGAAGAAGTCCGCGGCCAGGACACAGTTCCGGAAAGCACAGCCAGGATGCAGGGTGCCGGGAAAGCGCTGCATGAGTTGCTGCTGTCGGCGCAGCGTCAGGGCTGCCTCACTGCCGGCGTCTACGAGTCAGCCAAAGTCTTGAACGTGGACCCCGACAATGTGACCTTCTGTGTGCTGGCTGCGGGTGAGGAGGACGAGGGCGACATCGCGCTGCAGATCCATTTTACGCTGATCCAGGCTTTCTGCTGCGAGAACGACATCGACATAGTGCGCGTGGGCGATGTGCAGCGGCTGGCGGCTATCGTGGGCGCCGGCGAGGAGGCGGGTGCGCCGGGCGACCTGCACTGCATCCTCATTTCGAACCCCAACGAGGACGCCTGGAAGGATCCCGCCTTGGAGAAGCTCAGCCTGTTTTGCGAGGAGAGCCGCAGCGTTAACGACTGGGTGCCCAGCATCACCCTCCCCGAGTGA"

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
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    # The expected amino acid sequence for the HA tag
    ha_tag_aa = "YPYDVPDYA"
    
    # The final answer from the LLM to be checked
    llm_answer = "A" # The ribosome terminated the translation early

    protein_sequence = ""
    stop_codon_found = None
    stop_codon_position = -1

    # Ensure the sequence length is a multiple of 3 for clean processing
    # The final TGA is at the end, but we are looking for a *premature* one.
    seq_len = len(dna_sequence)
    
    # Translate the sequence
    for i in range(0, seq_len, 3):
        # Ensure we don't read past the end of the sequence
        if i + 3 > seq_len:
            break
        
        codon = dna_sequence[i:i+3]
        amino_acid = genetic_code.get(codon, '?') # '?' for unknown codons
        
        if amino_acid == '_':
            stop_codon_found = codon
            stop_codon_position = i
            break
        else:
            protein_sequence += amino_acid

    # --- Verification Steps ---

    # 1. Check for start codon and HA tag integrity (to evaluate option C)
    expected_start_and_tag = "M" + ha_tag_aa
    ha_tag_is_correct = protein_sequence.startswith(expected_start_and_tag)

    if not ha_tag_is_correct:
        return f"Incorrect. The answer is A, but the code found an issue with the HA tag itself. The expected start of the protein is '{expected_start_and_tag}', but the translated sequence begins with '{protein_sequence[:10]}'. This would support option C, not A."

    # 2. Check for premature termination (to evaluate option A)
    # The full DNA sequence is 360 bp, so a full protein would be ~120 aa.
    # A protein of length 11 is clearly truncated.
    is_premature_termination = stop_codon_found is not None and len(protein_sequence) < 20

    if not is_premature_termination:
        return f"Incorrect. The answer is A (early termination), but the code did not find a premature stop codon. The translated protein has a length of {len(protein_sequence)} amino acids."

    # 3. Consolidate findings
    # We found a premature stop codon, and the HA tag was correct.
    # The stop codon was TGA, which invalidates the premise of option B (mentions UAA/TAA).
    # The problem is at the translation level, making option D (proteolysis) a secondary effect at best, not the primary cause.
    # Therefore, the primary cause is indeed early termination.
    
    conclusion = "A" # Our code's conclusion matches option A

    if conclusion == llm_answer:
        return "Correct"
    else:
        return f"Incorrect. The provided answer was <<<{llm_answer}>>>, but the analysis points to a different conclusion. The code found a premature stop codon '{stop_codon_found}' at base position {stop_codon_position}, resulting in a truncated protein: '{protein_sequence}'. This strongly supports that the ribosome terminated translation early (Option A)."

# Run the check and print the result
print(check_answer())