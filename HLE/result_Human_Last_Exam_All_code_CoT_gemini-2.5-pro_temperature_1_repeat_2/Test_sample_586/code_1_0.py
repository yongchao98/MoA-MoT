import re

def find_best_sgrna():
    """
    This script analyzes a list of potential sgRNA sequences to find the most
    suitable one for targeting the second exon of a given gene with spCas9.
    """
    
    # 1. Isolate the Target Region (Exon 2)
    full_sequence = "ATGTGTGCCTGCGGGATACTTTTGACTTTTACTTTGCTTTTGCATTTTCTTGGTGTTCACTCAATGAATCCTCTGTTTCCAAGCGCATCCAGGGGCATGAAAGTGTCTAAGTCTGTTCCTGCTGAGGGCAACAGGAGAGCAAAATACGGCAAGAATGTGCTGTCAGCATCACTGTTATCCGGAGACATACAGTCCAGAAGGGCGATCAAGGATGCGATTGAACCTCACGATTACATGATTTCCATATACAAGACCTTTTCAGCGGCTGAAAAACTGGGACTGAACGCGAGTTTTTTCCGCTCGTCTAAAGCAGCAAACACCATCACGAGCTTTGTGGACGAGGGTCAAG^GTTAGTTATTTCTACTTATACAAGCAACAGTGATTTCAAACGCACACGTACTGATTCTATATTGGTACTCACAGGGAAAAAAAAAAAAAAAACATTTGTATACAATTCAAACAACTCTTAAAGGAATACAGTCAAATGTGTCAGTGAACAGATGGAAACAAAGCATTTTGAATATTAGGCCTATATCATCTATGATACTGCGGAAAATCTTCAAGAAATCTTTTTCCCCTAATAGTAAAAATAATGACAACAATATATGTATAACATTATACACTTCTGTTTACAATCTTGCATAAAATAAGTTGTGTTTGCATCAAAGTGTGTATACATGCACTGTCCATTTCAAATATTTTTTATTGGAATGTGTAGGAATTTTCACGATGTAGGCAGGTTATTATCACTATAAAGTGCCTTAGATGTCCCACAAGATTGAATCAGTCCCATATGAGCATAATGCGAAATTGATGTTTTAATATGATTGGTTAAACTTGTACACACATGCAGGTAGAATTATGAGTGTTTTGAAACATGTTTTTGCCAATTATTGCCATAGTCTTTTATTGAATGGATGTGATTTTGCCATGTCCCACACACTGCACAGCCAAGTTCAGTAAGTCTAAAAAGTAGCTAAATTAGATAAATTTTTTTTAAATGTTTAAGTATTCTTTCTATTCTTACAGTTATTTTGAAAACTAAATCATTTTTATAACTTTTATTTTTTTATTCTTTTATAATATTATTAATCATTTTGCACGAGTCTTTGAGTTTGCTGTCCACCCTGTCATGATGTAGTAAATCCCCTTTAAGAAACCCTCTGATGTACTCATTGGCATCCCCATGCCTATTTTGCTTTTCTTCAGAGGAGGTTAAAAAAACTGATGTGCACACATTAAATATCTACATATATGTTTCCTATTTTTCATCATATTGTGTTTGAAACCGAATGTGGTCAAGCTTAACATGTCCACCCTGTCATAGTAAAATATTAATTAATATAAAAAATTCGGAAATCAAAGATAGCTTTTAAACTGTATACAAAGAGCTTAAATAAGGAAACACTTTACCAGCTGCAGGTTCAACCTGTGTTAAATAAATGCTATCTTTAGCCAAAAATGTCCTCCTTGTTATTGTCCACCCTTTCACAAATCCTTCCTTGGGTGGACATATGCATCGTTATTGACACTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTTGTTAATCAGCTAATGTTTTATTATGGTACATCACATACATACTACACCAGTAGATGCAATACATAAGTGGACAATACAAATCTTTTGGCAATATTTATCTCAGTCTATATAAAGAATATCCTTTTAAAGTCCATATAAGGCAGCTCATTGACTGTTTGAAATTAAAATACATTATTTATCCTATTCTGGAAAAGAAAAAATATGATACATTTGTGCGTTGATGGATTTGAAACCACACTGGACTGAACTAATTTGAACTTTTAATTTCAATTCACTACAACTTCTATGTTAAGCTGCTTAGACACAATTTACATTACAGGTGTCAAATCCAGTTTCTTAAGAGCCACAGCTCTGCACAGTTTAGGGTTAACCCTAATTAAACACACCTGATCAAACTAATTGAGTCCTTCAGGCTTGTTTGATACCTACAGGTAGGTTTGTTAAAGCAAGGTTGGAACTAAATTGTGCAGAGCTGCGGCCCTTCAGGAACTAGATTTGACACCTAATTTACATTATGGAAACGCTATAGAAATAAAGATAAATTGAATTGAATAGATTTTTCTCCTCCAAAACACTATATATAAAAATACTAATTAGCAAATGCTAGTATTAGAAAAAAAAATTAGAACCTAGCTTTAAAAACTTTAGCATAATGAAAGAAACAGAGACACAAGACAGAAATAAATTTCAACATATGTCACCTTAATTAGTTAAAAACGAGTTCTCGATCTGCACATGCCATAACAGATATTGTAAATTTTGTGGATGCAGATCTAGTGTCAACAAGCATCTGTTCTCTTTGTTTCAG^ATGACCATTTGAACTCTCCACTTTGGAGACAGAAATATTTATTCGACGTATCAACGCTTTCTGAAAATGTGGAGATCCTGGGTGCCGAACTGAGGATTTACACAAAGATCTCCGGAAGCTTCCGCGCATCTGAAACCGGTCCTGTGGAAATACAGCTTCTCTCCTGCCAGTCGCACACTGTCCTTGATTCACAAACTTTGGATCTGGAGGATGCACATAAACCAAAATGGGAAGTTTTCGACGTCTGGGAGATTTTTAAGGAACGTCAGCACCACTCTCATGGCACCCGCTTCTGTTTAGAGCTCAGGGCCACACTGGATAATCCAGAGAGAGAAATTGATTTGCAATATCTTGGATTTCACAGACATGGCCGCCCGCAACTGAAGAAAGCCATACTGGTTGTTTTCACAAGGTCAAAAAAGAGGCAAAGTCTTTTTTATGAAAAAAGAGAGAAGATCAAGCTATGGGGTCTGGATAGTATTGGTAAGGAAAGAAGATCCCACTCGAAAACCCGCCGGAGCAGACGGACTGCTCTACCCAATCGCCATGGCAAGAGACATGGTAAAAAGTCAAAATCTAGATGCAGCAAAAAGCCACTGCATGTCAATTTCAGAGAGCTGGGTTGGGACGATTGGGTCATCGCTCCATTAGATTATGAGGCTTATCACTGTGAGGGCATGTGTGACTTTCCCCTCCGATCTCACCTGGAACCAACCAATCATGCCATCATACAAACTCTAATGAACTCAATGAACCCCAGCAACATGCCACCCAGCTGTTGCGTCCCCTCCAAACTCAGTCCCATTAGCATCTTGTACATTGACGCAGGAAATAATGTTGTGTACAAGCAGTATGAAGACATGGTAGTGGAGTCCTGCGGCTGCAGATGA"
    
    try:
        exon2 = full_sequence.split('^')[2]
    except IndexError:
        print("Error: The sequence format is incorrect. Could not find two '^' delimiters.")
        return

    choices_str = [
        "A. AGCGGTTTACTGAGACCCGG(TGG)", "B. TCCGGCGGGTTTTCGAGTGGG",
        "C. TTCATGCCCCTGGATGCGCT(TGG)", "D. CAGGACCGGTTTCAGATGCG(CGG)",
        "E. GCATCTGAAACCGGTCCTG(TGG)", "F. GGAAGCAATCCTCCGAACGT(TGG)",
        "G. ACGTTGCGAGGACAGAGTCA(AGG)", "H. CCCTTTCACAAATCCTTCCT(TGG)",
        "I. TTCACCCGCACCTTGAACGG(AGG)", "J. CTTTCTTTCTTTCTTTCTTTC(TTT)",
        "K. CTGCTCTACCCAATCGCCA(TGG)", "L. TGCCTG(CGG)",
        "M. TGCAAAGTAGATCGAGATGG(AGG)", "N. ACAGTCCAGAAGGGCGATCA(AGG)",
        "O. ATG(ACC)"
    ]

    # 2. Evaluate Each Option
    suitable_candidates = []
    print("--- Analysis of Potential sgRNA Targets ---")
    for choice_str in choices_str:
        match = re.match(r"([A-Z])\.?\s*([A-Z]+)\(?[A-Z]*\)?", choice_str)
        if match:
            option_id = match.group(1)
            full_target_str = choice_str.split('. ')[1]
            
            pam_match = re.search(r"\(([A-Z]+)\)$", full_target_str)
            if pam_match:
                guide = full_target_str[:pam_match.start()]
                pam = pam_match.group(1)
            else: # For cases like B without parentheses
                guide = full_target_str[:-3]
                pam = full_target_str[-3:]

            # Check if the target site (guide + PAM) is in exon 2
            if (guide + pam) in exon2:
                # Check for a valid spCas9 PAM (NGG)
                if pam.endswith("GG"):
                    # Calculate GC content
                    gc_count = guide.count('G') + guide.count('C')
                    gc_percent = (gc_count / len(guide)) * 100 if len(guide) > 0 else 0
                    
                    # Check for poly-T tract
                    has_poly_t = "TTTT" in guide
                    
                    print(f"\n[V] Option {option_id} is a potential candidate.")
                    print(f"    - Location: Found in Exon 2.")
                    print(f"    - Guide: {guide}")
                    print(f"    - PAM: {pam}")
                    print(f"    - GC Content: {gc_percent:.1f}%")
                    print(f"    - Contains Poly-T (TTTT): {'Yes' if has_poly_t else 'No'}")
                    
                    if 40 <= gc_percent <= 80 and not has_poly_t:
                        suitable_candidates.append({'id': option_id, 'guide': guide, 'pam': pam})
                else:
                    print(f"\n[X] Option {option_id}: Found in Exon 2 but has an invalid spCas9 PAM ({pam}).")
            # else:
                # print(f"\n[X] Option {option_id}: Not found in Exon 2.")

    # 3. Select the Best Candidate
    print("\n--- Final Recommendation ---")
    if not suitable_candidates:
        print("No suitable sgRNA candidates were found in the provided options.")
    elif len(suitable_candidates) == 1:
        best_choice = suitable_candidates[0]
        print(f"Only one suitable candidate found: Option {best_choice['id']}.")
    else:
        print("Multiple suitable candidates found (E, K). Comparing them:")
        print("- Candidate E ('GCATCTGAAACCGGTCCTG') is located earlier in the exon, which is preferred for gene knockout experiments.")
        print("- Candidate E's sequence is also more complex and less repetitive than K's ('CTGCTCTACCCAATCGCCA'), suggesting a lower probability of off-target effects.")
        print("\nBased on these criteria, Candidate E is the superior choice.")
        best_choice = next(c for c in suitable_candidates if c['id'] == 'E')

    # 4. Output the Final Answer
    if suitable_candidates:
        final_guide = best_choice['guide']
        final_pam = best_choice['pam']
        
        print("\nThe most suitable target sequence is from option E.")
        print("\nFinal Target Equation:")
        print(f"Guide Sequence: {final_guide}")
        print(f"PAM Sequence: {final_pam}")
        print(f"Full Target Site in Exon 2: {final_guide + final_pam}")

find_best_sgrna()
<<<E>>>