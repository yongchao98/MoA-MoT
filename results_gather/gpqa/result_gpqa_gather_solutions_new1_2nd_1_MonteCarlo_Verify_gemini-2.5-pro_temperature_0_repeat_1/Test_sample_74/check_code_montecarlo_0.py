import sys
from io import StringIO

def check_answer():
    """
    This function checks the correctness of the provided answer by analyzing the DNA sequence.
    It simulates the process of translation to identify the molecular reason for the failed protein expression.
    """
    
    # Store the original stdout to restore it later
    original_stdout = sys.stdout
    # Redirect stdout to a string buffer to capture prints
    captured_output = StringIO()
    sys.stdout = captured_output

    # --- Step 1: Define the problem's parameters ---
    dna_sequence = "ATGTACCCATACGATGTTCCAGATTACGCCAAATGACTCTGGAAGAAGTCCGCGGCCAGGACACAGTTCCGGAAAGCACAGCCAGGATGCAGGGTGCCGGGAAAGCGCTGCATGAGTTGCTGCTGTCGGCGCAGCGTCAGGGCTGCCTCACTGCCGGCGTCTACGAGTCAGCCAAAGTCTTGAACGTGGACCCCGACAATGTGACCTTCTGTGTGCTGGCTGCGGGTGAGGAGGACGAGGGCGACATCGCGCTGCAGATCCATTTTACGCTGATCCAGGCTTTCTGCTGCGAGAACGACATCGACATAGTGCGCGTGGGCGATGTGCAGCGGCTGGCGGCTATCGTGGGCGCCGGCGAGGAGGCGGGTGCGCCGGGCGACCTGCACTGCATCCTCATTTCGAACCCCAACGAGGACGCCTGGAAGGATCCCGCCTTGGAGAAGCTCAGCCTGTTTTGCGAGGAGAGCCGCAGCGTTAACGACTGGGTGCCCAGCATCACCCTCCCCGAGTGA"
    
    # The options as presented in the question
    options = {
        "A": "The tRNA for the UAA codon does not exist in the mouse",
        "B": "The sequence for the antigenic determinant has a missense mutation",
        "C": "The ribosome terminated the translation early",
        "D": "The lack of the linker sequence is triggering proteolysis of the nascent chain"
    }
    
    # The provided answer to check
    llm_answer = "C"

    # Standard DNA codon table
    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'AGA':'R', 'AGG':'R', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'TTA':'L', 'TTG':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GGC':'G', 'GGA':'G', 'GGG':'G', 'GGT':'G',
        'TTC':'F', 'TTT':'F', 'TGG':'W', 'TAC':'Y', 'TAT':'Y',
        'GTC':'V', 'GTA':'V', 'GTG':'V', 'GTT':'V',
        'TGC':'C', 'TGT':'C', 'TAA':'Stop', 'TAG':'Stop', 'TGA':'Stop'
    }

    print("--- Verification Process ---")

    # --- Step 2: Translate the DNA sequence ---
    protein_sequence = ""
    stop_codon_found = None
    stop_codon_position = -1
    
    # Check if the sequence starts with a start codon
    if dna_sequence.startswith("ATG"):
        print("Analysis: Sequence starts with a valid start codon 'ATG'.")
        # Start translation from the first codon
        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i+3]
            if len(codon) == 3:
                amino_acid = genetic_code.get(codon, '?')
                if amino_acid == 'Stop':
                    stop_codon_found = codon
                    stop_codon_position = i
                    print(f"Analysis: Found a '{codon}' stop codon at base position {i+1}.")
                    break
                protein_sequence += amino_acid
    else:
        print("Analysis: Sequence does not start with a valid 'ATG' start codon.")


    # --- Step 3: Evaluate each option based on the analysis ---
    
    # Check Option A: tRNA for UAA
    print("\n--- Evaluating Option A ---")
    print(f"Claim: {options['A']}")
    print("Analysis: This claim is biologically incorrect. Stop codons (like UAA, UAG, UGA) are recognized by protein release factors, not tRNAs.")
    print(f"Analysis: Furthermore, the stop codon found in the sequence was '{stop_codon_found}', not 'TAA' (which corresponds to mRNA's UAA).")
    is_A_correct = False
    print("Conclusion: Option A is incorrect.")

    # Check Option B: Missense mutation in HA tag
    print("\n--- Evaluating Option B ---")
    print(f"Claim: {options['B']}")
    expected_ha_tag_aa = "YPYDVPDYA"
    # The HA tag DNA is the 27 bases (9 codons) following the initial ATG start codon.
    ha_tag_dna = dna_sequence[3:3+27]
    translated_ha_tag = ""
    for i in range(0, len(ha_tag_dna), 3):
        codon = ha_tag_dna[i:i+3]
        translated_ha_tag += genetic_code.get(codon, '?')
    
    print(f"Analysis: Expected HA tag amino acid sequence: {expected_ha_tag_aa}")
    print(f"Analysis: Translated HA tag from DNA sequence: {translated_ha_tag}")
    if translated_ha_tag == expected_ha_tag_aa:
        print("Conclusion: Option B is incorrect. The HA tag sequence is correct.")
        is_B_correct = False
    else:
        print("Conclusion: Option B is correct. There is a mutation in the HA tag.")
        is_B_correct = True

    # Check Option C: Ribosome terminated early
    print("\n--- Evaluating Option C ---")
    print(f"Claim: {options['C']}")
    expected_full_protein_len = len(dna_sequence) // 3
    actual_protein_len = len(protein_sequence)
    print(f"Analysis: A stop codon was found at base position {stop_codon_position+1}.")
    print(f"Analysis: This resulted in a truncated peptide of only {actual_protein_len} amino acids.")
    print(f"Analysis: The full gene could have coded for up to {expected_full_protein_len} amino acids.")
    if stop_codon_position != -1 and actual_protein_len < 50: # Using < 50 as a heuristic for "early"
        print("Conclusion: Option C is correct. The presence of a premature stop codon caused early termination.")
        is_C_correct = True
    else:
        print("Conclusion: Option C is incorrect.")
        is_C_correct = False

    # Check Option D: Proteolysis due to lack of linker
    print("\n--- Evaluating Option D ---")
    print(f"Claim: {options['D']}")
    print("Analysis: Proteolysis is the degradation of a protein *after* it has been synthesized.")
    if is_C_correct:
        print("Analysis: Our findings show that the full-length protein was never synthesized due to a premature stop codon.")
        print("Conclusion: Option D is incorrect as the primary issue is a failure of synthesis, not degradation.")
        is_D_correct = False
    else:
        print("Conclusion: This could be a possibility, but a more direct cause was not found.")
        is_D_correct = False # Cannot be confirmed by sequence analysis alone, but is secondary to synthesis failure.

    # --- Step 4: Final Verdict ---
    print("\n--- Final Verdict ---")
    correct_option = None
    if is_A_correct: correct_option = "A"
    elif is_B_correct: correct_option = "B"
    elif is_C_correct: correct_option = "C"
    elif is_D_correct: correct_option = "D"

    print(f"Programmatic analysis determined the correct option is '{correct_option}'.")
    print(f"The LLM's answer was '{llm_answer}'.")

    # Restore stdout
    sys.stdout = original_stdout
    # Get the captured output as a string
    result_string = captured_output.getvalue()

    if llm_answer == correct_option:
        return "Correct"
    else:
        return f"Incorrect. The provided answer is {llm_answer}, but the analysis shows {correct_option} is the correct choice.\n\n--- Detailed Analysis Log ---\n{result_string}"

# Run the check and print the result
result = check_answer()
print(result)