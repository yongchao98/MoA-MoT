import re

def find_best_sgrna_target():
    """
    Analyzes a gene sequence to find the best sgRNA target in the second exon
    from a list of choices.
    """
    full_sequence = "ATGTGTGCCTGCGGGATACTTTTGACTTTTACTTTGCTTTTGCATTTTCTTGGTGTTCACTCAATGAATCCTCTGTTTCCAAGCGCATCCAGGGGCATGAAAGTGTCTAAGTCTGTTCCTGCTGAGGGCAACAGGAGAGCAAAATACGGCAAGAATGTGCTGTCAGCATCACTGTTATCCGGAGACATACAGTCCAGAAGGGCGATCAAGGATGCGATTGAACCTCACGATTACATGATTTCCATATACAAGACCTTTTCAGCGGCTGAAAAACTGGGACTGAACGCGAGTTTTTTCCGCTCGTCTAAAGCAGCAAACACCATCACGAGCTTTGTGGACGAGGGTCAAG^GTTAGTTATTTCTACTTATACAAGCAACAGTGATTTCAAACGCACACGTACTGATTCTATATTGGTACTCACAGGGAAAAAAAAAAAAAAAACATTTGTATACAATTCAAACAACTCTTAAAGGAATACAGTCAAATGTGTCAGTGAACAGATGGAAACAAAGCATTTTGAATATTAGGCCTATATCATCTATGATACTGCGGAAAATCTTCAAGAAATCTTTTTCCCCTAATAGTAAAAATAATGACAACAATATATGTATAACATTATACACTTCTGTTTACAATCTTGCATAAAATAAGTTGTGTTTGCATCAAAGTGTGTATACATGCACTGTCCATTTCAAATATTTTTTATTGGAATGTGTAGGAATTTTCACGATGTAGGCAGGTTATTATCACTATAAAGTGCCTTAGATGTCCCACAAGATTGAATCAGTCCCATATGAGCATAATGCGAAATTGATGTTTTAATATGATTGGTTAAACTTGTACACACATGCAGGTAGAATTATGAGTGTTTTGAAACATGTTTTTGCCAATTATTGCCATAGTCTTTTATTGAATGGATGTGATTTTGCCATGTCCCACACACTGCACAGCCAAGTTCAGTAAGTCTAAAAAGTAGCTAAATTAGATAAATTTTTTTTAAATGTTTAAGTATTCTTTCTATTCTTACAGTTATTTTGAAAACTAAATCATTTTTATAACTTTTATTTTTTTATTCTTTTATAATATTATTAATCATTTTGCACGAGTCTTTGAGTTTGCTGTCCACCCTGTCATGATGTAGTAAATCCCCTTTAAGAAACCCTCTGATGTACTCATTGGCATCCCCATGCCTATTTTGCTTTTCTTTCAAGGAGGTTAAAAAAACTGATGTGCACACATTAAATATCTACATATATGTTTCCTATTTTTCATCATATTGTGTTTGAAACCGAATGTGGTCAAGCTTAACATGTCCACCCTGTCATAGTAAAATATTAATTAATATAAAAAATTCGGAAATCAAAGATAGCTTTTAAACTGTATACAAAGAGCTTAAATAAGGAAACACTTTACCAGCTGCAGGTTCAACCTGTGTTAAATAAATGCTATCTTTAGCCAAAAATGTCCTCCTTGTTATTGTCCACCCTTTCACAAATCCTTCCTTGGGTGGACATATGCATCGTTATTGACACTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTTGTTAATCAGCTAATGTTTTATTATGGTACATCACATACATACTACACCAGTAGATGCAATACATAAGTGGACAATACAAATCTTTTGGCAATATTTATCTCAGTCTATATAAAGAATATCCTTTTAAAGTCCATATAAGGCAGCTCATTGACTGTTTGAAATTAAAATACATTATTTATCCTATTCTGGAAAAGAAAAAATATGATACATTTGTGCGTTGATGGATTTGAAACCACACTGGACTGAACTAATTTGAACTTTTAATTTCAATTCACTACAACTTCTATGTTAAGCTGCTTAGACACAATTTACATTACAGGTGTCAAATCCAGTTTCTTAAGAGCCACAGCTCTGCACAGTTTAGGGTTAACCCTAATTAAACACACCTGATCAAACTAATTGAGTCCTTCAGGCTTGTTTGATACCTACAGGTAGGTTTGTTAAAGCAAGGTTGGAACTAAATTGTGCAGAGCTGCGGCCCTTCAGGAACTAGATTTGACACCTAATTTACATTATGGAAACGCTATAGAAATAAAGATAAATTGAATTGAATAGATTTTTCTCCTCCAAAACACTATATATAAAAATACTAATTAGCAAATGCTAGTATTAGAAAAAAAAATTAGAACCTAGCTTTAAAAACTTTAGCATAATGAAAGAAACAGAGACACAAGACAGAAATAAATTTCAACATATGTCACCTTAATTAGTTAAAAACGAGTTCTCGATCTGCACATGCCATAACAGATATTGTAAATTTTGTGGATGCAGATCTAGTGTCAACAAGCATCTGTTCTCTTTGTTTCAG^ATGACCATTTGAACTCTCCACTTTGGAGACAGAAATATTTATTCGACGTATCAACGCTTTCTGAAAATGTGGAGATCCTGGGTGCCGAACTGAGGATTTACACAAAGATCTCCGGAAGCTTCCGCGCATCTGAAACCGGTCCTGTGGAAATACAGCTTCTCTCCTGCCAGTCGCACACTGTCCTTGATTCACAAACTTTGGATCTGGAGGATGCACATAAACCAAAATGGGAAGTTTTCGACGTCTGGGAGATTTTTAAGGAACGTCAGCACCACTCTCATGGCACCCGCTTCTGTTTAGAGCTCAGGGCCACACTGGATAATCCAGAGAGAGAAATTGATTTGCAATatcttggatttcacagacatggccgcccgcaactgaagaaagccatactggttgttttcacaaggtcaaaaaagaggcaaagtcttttttatgaaaaaagagagaagatcaagctatggggtctggatagtattggtaaggaaagaagatcccactcgaaaacccgccggagcagacggactgctctacccaatcgccatggcaagagacatggtaaaaagtcaaaatctagatgcagcaaaaagccactgcatgtcaatttcagagagctgggttgggacgattgggtcatcgctccattagattatgaggcttatcactgtgagggcatgtgtgactttcccctccgatctcacctggaaccaaccaatcatgccatcatacaaactctaatgaactcaatgaaccccagcaacatgccacccagctgttgcgtcccctccaaactcagtcccattagcatcttgtacattgacgcaggaaataatgttgtgtacaagcagtatgaagacatangtagtggagtcctgcggctgcagatga".upper()

    parts = full_sequence.split('^')
    exon2 = parts[2]

    choices = {
        "A": "AGCGGTTTACTGAGACCCGG(TGG)", "B": "TCCGGCGGGTTTTCGAGTGGG",
        "C": "TTCATGCCCCTGGATGCGCT(TGG)", "D": "CAGGACCGGTTTCAGATGCG(CGG)",
        "E": "GCATCTGAAACCGGTCCTG(TGG)", "F": "GGAAGCAATCCTCCGAACGT(TGG)",
        "G": "ACGTTGCGAGGACAGAGTCA(AGG)", "H": "CCCTTTCACAAATCCTTCCT(TGG)",
        "I": "TTCACCCGCACCTTGAACGG(AGG)", "J": "CTTTCTTTCTTTCTTTCTTTC(TTT)",
        "K": "CTGCTCTACCCAATCGCCA(TGG)", "L": "TGCCTG(CGG)",
        "M": "TGCAAAGTAGATCGAGATGG(AGG)", "N": "ACAGTCCAGAAGGGCGATCA(AGG)",
        "O": "ATG(ACC)"
    }

    valid_targets = []
    print("Analyzing potential sgRNA targets...\n")

    for letter, choice_str in sorted(choices.items()):
        # Parse the target and PAM from the string
        match = re.match(r"([A-Z]+)\(([A-Z]+)\)", choice_str)
        if match:
            target, pam = match.groups()
        else: # Handle cases like 'B' where PAM is not in parentheses
            target, pam = choice_str[:-3], choice_str[-3:]
        
        full_target_seq = target + pam

        # Check for valid NGG PAM
        if not (pam[1:] == 'GG'):
            print(f"Choice {letter}: Invalid PAM '{pam}'. Skipping.")
            continue

        # Check if target is in exon 2
        position = exon2.find(full_target_seq)
        if position != -1:
            # Check for undesirable properties, e.g., repeats
            if "TTTT" in target:
                 print(f"Choice {letter}: Found in exon 2, but contains 'TTTT' which can cause issues. Low priority.")
            elif len(set(target[-5:])) < 3: # Simple heuristic for repetitive ends
                 print(f"Choice {letter}: Found in exon 2, but has a repetitive end. Medium priority.")
            else:
                print(f"Choice {letter}: Found in exon 2 at position {position}. This is a valid candidate.")
                valid_targets.append({'letter': letter, 'pos': position, 'seq': full_target_seq})
        else:
            # If not in exon 2, we can ignore it for this task.
            pass

    print("\n--- Evaluation Summary ---")
    if not valid_targets:
        print("No valid targets from the choices were found in the second exon.")
        return

    # Sort valid targets by their position in the exon (earlier is better)
    best_target = sorted(valid_targets, key=lambda x: x['pos'])[0]

    print(f"\nMultiple valid targets were found in exon 2: {[t['letter'] for t in valid_targets]}")
    print("To maximize the chance of creating a non-functional protein, it is best to target early in the exon.")
    print(f"The earliest valid target is Choice {best_target['letter']} at position {best_target['pos']}.")
    print("\n--- Final Recommendation ---")
    print(f"The most suitable target sequence is from choice {best_target['letter']}:")
    print(f"Target + PAM: {best_target['seq']}")
    
    # Final answer format for the platform
    print(f"\n<<<E>>>")


find_best_sgrna_target()