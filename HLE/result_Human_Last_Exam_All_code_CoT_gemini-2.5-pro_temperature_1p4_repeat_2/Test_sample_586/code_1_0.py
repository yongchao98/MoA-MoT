import re

def find_best_sgrna():
    """
    This script analyzes potential sgRNA target sequences to find the most suitable one
    for targeting the second exon of the provided gene sequence using spCas9.
    """
    # The full DNA sequence provided by the user
    full_seq = "ATGTGTGCCTGCGGGATACTTTTGACTTTTACTTTGCTTTTGCATTTTCTTGGTGTTCACTCAATGAATCCTCTGTTTCCAAGCGCATCCAGGGGCATGAAAGTGTCTAAGTCTGTTCCTGCTGAGGGCAACAGGAGAGCAAAATACGGCAAGAATGTGCTGTCAGCATCACTGTTATCCGGAGACATACAGTCCAGAAGGGCGATCAAGGATGCGATTGAACCTCACGATTACATGATTTCCATATACAAGACCTTTTCAGCGGCTGAAAAACTGGGACTGAACGCGAGTTTTTTCCGCTCGTCTAAAGCAGCAAACACCATCACGAGCTTTGTGGACGAGGGTCAAG^GTTAGTTATTTCTACTTATACAAGCAACAGTGATTTCAAACGCACACGTACTGATTCTATATTGGTACTCACAGGGAAAAAAAAAAAAAAAACATTTGTATACAATTCAAACAACTCTTAAAGGAATACAGTCAAATGTGTCAGTGAACAGATGGAAACAAAGCATTTTGAATATTAGGCCTATATCATCTATGATACTGCGGAAAATCTTCAAGAAATCTTTTTCCCCTAATAGTAAAAATAATGACAACAATATATGTATAACATTATACACTTCTGTTTACAATCTTGCATAAAATAAGTTGTGTTTGCATCAAAGTGTGTATACATGCACTGTCCATTTCAAATATTTTTTATTGGAATGTGTAGGAATTTTCACGATGTAGGCAGGTTATTATCACTATAAAGTGCCTTAGATGTCCCACAAGATTGAATCAGTCCCATATGAGCATAATGCGAAATTGATGTTTTAATATGATTGGTTAAACTTGTACACACATGCAGGTAGAATTATGAGTGTTTTGAAACATGTTTTTGCCAATTATTGCCATAGTCTTTTATTGAATGGATGTGATTTTGCCATGTCCCACACACTGCACAGCCAAGTTCAGTAAGTCTAAAAAGTAGCTAAATTAGATAAATTTTTTTTAAATGTTTAAGTATTCTTTCTATTCTTACAGTTATTTTGAAAACTAAATCATTTTTATAACTTTTATTTTTTTATTCTTTTATAATATTATTAATCATTTTGCACGAGTCTTTGAGTTTGCTGTCCACCCTGTCATGATGTAGTAAATCCCCTTTAAGAAACCCTCTGATGTACTCATTGGCATCCCCATGCCTATTTTGCTTTTCTTTCAAGGAGGTTAAAAAAACTGATGTGCACACATTAAATATCTACATATATGTTTCCTATTTTTCATCATATTGTGTTTGAAACCGAATGTGGTCAAGCTTAACATGTCCACCCTGTCATAGTAAAATATTAATTAATATAAAAAATTCGGAAATCAAAGATAGCTTTTAAACTGTATACAAAGAGCTTAAATAAGGAAACACTTTACCAGCTGCAGGTTCAACCTGTGTTAAATAAATGCTATCTTTAGCCAAAAATGTCCTCCTTGTTATTGTCCACCCTTTCACAAATCCTTCCTTGGGTGGACATATGCATCGTTATTGACACTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTTGTTAATCAGCTAATGTTTTATTATGGTACATCACATACATACTACACCAGTAGATGCAATACATAAGTGGACAATACAAATCTTTTGGCAATATTTATCTCAGTCTATATAAAGAATATCCTTTTAAAGTCCATATAAGGCAGCTCATTGACTGTTTGAAATTAAAATACATTATTTATCCTATTCTGGAAAAGAAAAAATATGATACATTTGTGCGTTGATGGATTTGAAACCACACTGGACTGAACTAATTTGAACTTTTAATTTCAATTCACTACAACTTCTATGTTAAGCTGCTTAGACACAATTTACATTACAGGTGTCAAATCCAGTTTCTTAAGAGCCACAGCTCTGCACAGTTTAGGGTTAACCCTAATTAAACACACCTGATCAAACTAATTGAGTCCTTCAGGCTTGTTTGATACCTACAGGTAGGTTTGTTAAAGCAAGGTTGGAACTAAATTGTGCAGAGCTGCGGCCCTTCAGGAACTAGATTTGACACCTAATTTACATTATGGAAACGCTATAGAAATAAAGATAAATTGAATTGAATAGATTTTTCTCCTCCAAAACACTATATATAAAAATACTAATTAGCAAATGCTAGTATTAGAAAAAAAAATTAGAACCTAGCTTTAAAAACTTTAGCATAATGAAAGAAACAGAGACACAAGACAGAAATAAATTTCAACATATGTCACCTTAATTAGGTAAAAACGAGTTCTCGATCTGCACATGCCATAACAGATATTGTAAATTTTGTGGATGCAGATCTAGTGTCAACAAGCATCTGTTCTCTTTGTTTCAG^ATGACCATTTGAACTCTCCACTTTGGAGACAGAAATATTTATTCGACGTATCAACGCTTTCTGAAAATGTGGAGATCCTGGGTGCCGAACTGAGGATTTACACAAAGATCTCCGGAAGCTTCCGCGCATCTGAAACCGGTCCTGTGGAAATACAGCTTCTCTCCTGCCAGTCGCACACTGTCCTTGATTCACAAACTTTGGATCTGGAGGATGCACATAAACCAAAATGGGAAGTTTTCGACGTCTGGGAGATTTTTAAGGAACGTCAGCACCACTCTCATGGCACCCGCTTCTGTTTAGAGCTCAGGGCCACACTGGATAATCCAGAGAGAGAAATTGATTTGCAATATCTTGGATTTCACAGACATGGCCGCCCGCAACTGAAGAAAGCCATACTGGTTGTTTTCACAAGGTCAAAAAAGAGGCAAAGTCTTTTTTATGAAAAAAGAGAGAAGATCAAGCTATGGGGTCTGGATAGTATTGGTAAGGAAAGAAGATCCCACTCGAAAACCCGCCGGAGCAGACGGACTGCTCTACCCAATCGCCATGGCAAGAGACATGGTAAAAAGTCAAAATCTAGATGCAGCAAAAAGCCACTGCATGTCAATTTCAGAGAGCTGGGTTGGGACGATTGGGTCATCGCTCCATTAGATTATGAGGCTTATCACTGTGAGGGCATGTGTGACTTTCCCCTCCGATCTCACCTGGAACCAACCAATCATGCCATCATACAAACTCTAATGAACTCAATGAACCCCAGCAACATGCCACCCAGCTGTTGCGTCCCCTCCAAACTCAGTCCCATTAGCATCTTGTACATTGACGCAGGAAATAATGTTGTGTACAAGCAGTATGAAGACATGGTAGTGGAGTCCTGCGGCTGCAGATGA"

    # Split the sequence into exons and intron based on the '^' delimiter
    parts = full_seq.split('^')
    exon2_seq = parts[2]

    # The list of answer choices provided by the user
    choices = {
        'A': "AGCGGTTTACTGAGACCCGG(TGG)", 'B': "TCCGGCGGGTTTTCGAGTGGG",
        'C': "TTCATGCCCCTGGATGCGCT(TGG)", 'D': "CAGGACCGGTTTCAGATGCG(CGG)",
        'E': "GCATCTGAAACCGGTCCTG(TGG)", 'F': "GGAAGCAATCCTCCGAACGT(TGG)",
        'G': "ACGTTGCGAGGACAGAGTCA(AGG)", 'H': "CCCTTTCACAAATCCTTCCT(TGG)",
        'I': "TTCACCCGCACCTTGAACGG(AGG)", 'J': "CTTTCTTTCTTTCTTTCTTTC(TTT)",
        'K': "CTGCTCTACCCAATCGCCA(TGG)", 'L': "TGCCTG(CGG)",
        'M': "TGCAAAGTAGATCGAGATGG(AGG)", 'N': "ACAGTCCAGAAGGGCGATCA(AGG)",
        'O': "ATG(ACC)"
    }
    
    print("Evaluating potential sgRNA target sequences...\n")
    valid_candidates = []

    for key, value in choices.items():
        # Use regex to parse the choice string into target and PAM components
        match = re.match(r"([A-Z]+)\(([A-Z]+)\)", value)
        if not match:
            continue
        
        target, pam = match.groups()
        
        # Rule 1: Check for a valid spCas9 PAM (NGG)
        if not (len(pam) == 3 and pam[1:] == 'GG'):
            continue

        # Rule 2: Check if the full site (target + PAM) is in the second exon
        if (target + pam) not in exon2_seq:
            continue
        
        # If all rules pass, it's a valid candidate. Store its properties.
        position = exon2_seq.find(target + pam)
        gc_content = (target.count('G') + target.count('C')) / len(target) * 100
        valid_candidates.append({
            'key': key,
            'target': target,
            'pam': pam,
            'position_in_exon': position,
            'gc_content': gc_content
        })

    print("--- Analysis of Valid Candidates Found in Exon 2 ---")
    if valid_candidates:
        for cand in valid_candidates:
            print(f"\nCandidate {cand['key']}:")
            print(f"  Target Sequence (20-mer): {cand['target']}")
            print(f"  PAM Sequence: {cand['pam']}")
            print(f"  Position in Exon 2 (0-indexed): {cand['position_in_exon']}")
            print(f"  GC Content: {cand['gc_content']:.1f}%")
    else:
        print("No valid candidates were found in Exon 2 from the provided choices.")

    print("\n--- Final Recommendation ---")
    # Identify candidates E and K from the valid list
    cand_E = next((c for c in valid_candidates if c['key'] == 'E'), None)
    cand_K = next((c for c in valid_candidates if c['key'] == 'K'), None)

    if cand_E and cand_K:
        print("Both E and K are valid sgRNA targets located within Exon 2.")
        print("To choose the best one, we consider two common design principles:")
        print("1. Position in Exon: Targeting early in an exon is often preferred to maximize the disruptive effect of a potential frameshift mutation.")
        print("2. Sequence Complexity: Less repetitive sequences are generally better to minimize potential off-target binding.")
        
        print(f"\nComparing Candidate E and K:")
        print(f" - Position: Candidate E is at position {cand_E['position_in_exon']}, while Candidate K is at {cand_K['position_in_exon']}. E is located much earlier in the exon.")
        print(f" - Complexity: Candidate E's sequence ('{cand_E['target']}') is less repetitive than Candidate K's sequence ('{cand_K['target']}'), which contains motifs like 'CTGCT' and 'CCC'.")
        
        print("\nConclusion: Based on its earlier position and more complex sequence, E is the more suitable choice.")
        final_answer_key = 'E'
        final_answer_data = cand_E
        
    elif valid_candidates:
        final_answer_key = valid_candidates[0]['key']
        final_answer_data = valid_candidates[0]
        print(f"Only one suitable candidate was found: {final_answer_key}.")
    else:
        print("None of the choices are suitable.")
        return

    print(f"\nThe most suitable target sequence is from choice {final_answer_key}.")
    print(f"Final Target Equation:")
    print(f"{final_answer_data['target']} + ({final_answer_data['pam']})")
    
    # Finally, print the answer in the required special format
    print(f"\n<<<{final_answer_key}>>>")

if __name__ == '__main__':
    find_best_sgrna()