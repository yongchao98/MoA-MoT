import re

def analyze_sgrna_targets():
    """
    Analyzes a list of potential sgRNA targets against a given DNA sequence,
    focusing on the second exon.
    """
    # The full DNA sequence provided by the user
    full_sequence = "ATGTGTGCCTGCGGGATACTTTTGACTTTTACTTTGCTTTTGCATTTTCTTGGTGTTCACTCAATGAATCCTCTGTTTCCAAGCGCATCCAGGGGCATGAAAGTGTCTAAGTCTGTTCCTGCTGAGGGCAACAGGAGAGCAAAATACGGCAAGAATGTGCTGTCAGCATCACTGTTATCCGGAGACATACAGTCCAGAAGGGCGATCAAGGATGCGATTGAACCTCACGATTACATGATTTCCATATACAAGACCTTTTCAGCGGCTGAAAAACTGGGACTGAACGCGAGTTTTTTCCGCTCGTCTAAAGCAGCAAACACCATCACGAGCTTTGTGGACGAGGGTCAAG^GTTAGTTATTTCTACTTATACAAGCAACAGTGATTTCAAACGCACACGTACTGATTCTATATTGGTACTCACAGGGAAAAAAAAAAAAAAAACATTTGTATACAATTCAAACAACTCTTAAAGGAATACAGTCAAATGTGTCAGTGAACAGATGGAAACAAAGCATTTTGAATATTAGGCCTATATCATCTATGATACTGCGGAAAATCTTCAAGAAATCTTTTTCCCCTAATAGTAAAAATAATGACAACAATATATGTATAACATTATACACTTCTGTTTACAATCTTGCATAAAATAAGTTGTGTTTGCATCAAAGTGTGTATACATGCACTGTCCATTTCAAATATTTTTTATTGGAATGTGTAGGAATTTTCACGATGTAGGCAGGTTATTATCACTATAAAGTGCCTTAGATGTCCCACAAGATTGAATCAGTCCCATATGAGCATAATGCGAAATTGATGTTTTAATATGATTGGTTAAACTTGTACACACATGCAGGTAGAATTATGAGTGTTTTGAAACATGTTTTTGCCAATTATTGCCATAGTCTTTTATTGAATGGATGTGATTTTGCCATGTCCCACACACTGCACAGCCAAGTTCAGTAAGTCTAAAAAGTAGCTAAATTAGATAAATTTTTTTTAAATGTTTAAGTATTCTTTCTATTCTTACAGTTATTTTGAAAACTAAATCATTTTTATAACTTTTATTTTTTTATTCTTTTATAATATTATTAATCATTTTGCACGAGTCTTTGAGTTTGCTGTCCACCCTGTCATGATGTAGTAAATCCCCTTTAAGAAACCCTCTGATGTACTCATTGGCATCCCCATGCCTATTTTGCTTTTCTTCAGAGGAGGTTAAAAAAACTGATGTGCACACATTAAATATCTACATATATGTTTCCTATTTTTCATCATATTGTGTTTGAAACCGAATGTGGTCAAGCTTAACATGTCCACCCTGTCATAGTAAAATATTAATTAATATAAAAAATTCGGAAATCAAAGATAGCTTTTAAACTGTATACAAAGAGCTTAAATAAGGAAACACTTTACCAGCTGCAGGTTCAACCTGTGTTAAATAAATGCTATCTTTAGCCAAAAATGTCCTCCTTGTTATTGTCCACCCTTTCACAAATCCTTCCTTGGGTGGACATATGCATCGTTATTGACACTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTTGTTAATCAGCTAATGTTTTATTATGGTACATCACATACATACTACACCAGTAGATGCAATACATAAGTGGACAATACAAATCTTTTGGCAATATTTATCTCAGTCTATATAAAGAATATCCTTTTAAAGTCCATATAAGGCAGCTCATTGACTGTTTGAAATTAAAATACATTATTTATCCTATTCTGGAAAAGAAAAAATATGATACATTTGTGCGTTGATGGATTTGAAACCACACTGGACTGAACTAATTTGAACTTTTAATTTCAATTCACTACAACTTCTATGTTAAGCTGCTTAGACACAATTTACATTACAGGTGTCAAATCCAGTTTCTTAAGAGCCACAGCTCTGCACAGTTTAGGGTTAACCCTAATTAAACACACCTGATCAAACTAATTGAGTCCTTCAGGCTTGTTTGATACCTACAGGTAGGTTTGTTAAAGCAAGGTTGGAACTAAATTGTGCAGAGCTGCGGCCCTTCAGGAACTAGATTTGACACCTAATTTACATTATGGAAACGCTATAGAAATAAAGATAAATTGAATTGAATAGATTTTTCTCCTCCAAAACACTATATATAAAAATACTAATTAGCAAATGCTAGTATTAGAAAAAAAAATTAGAACCTAGCTTTAAAAACTTTAGGATAATGAAAGAAACAGAGACACAAGACAGAAATAAATTTCAACATATGTCACCTTAATTAGTTAAAAACGAGTTCTCGATCTGCACATGCCATAACAGATATTGTAAattttgtggatgcagatctagtgtcaacaagcatctgttctctttgtttcag^ATGACCATTTGAACTCTCCACTTTGGAGACAGAAATATTTATTCGACGTATCAACGCTTTCTGAAAATGTGGAGATCCTGGGTGCCGAACTGAGGATTTACACAAAGATCTCCGGAAGCTTCCGCGCATCTGAAACCGGTCCTGTGGAAATACAGCTTCTCTCCTGCCAGTCGCACACTGTCCTTGATTCACAAACTTTGGATCTGGAGGATGCACATAAACCAAAATGGGAAGTTTTCGACGTCTGGGAGATTTTTAAGGAACGTCAGCACCACTCTCATGGCACCCGCTTCTGTTTAGAGCTCAGGGCCACACTGGATAATCCAGAGAGAGAAATTGATTTGCAATATCTTGGATTTCACAGACATGGCCGCCCGCAACTGAAGAAAGCCATACTGGTTGTTTTCACAAGGTCAAAAAAGAGGCAAAGTCTTTTTTATGAAAAAAGAGAGAAGATCAAGCTATGGGGTCTGGATAGTATTGGTAAGGAAAGAAGATCCCACTCGAAAACCCGCCGGAGCAGACGGACTGCTCTACCCAATCGCCATGGCAAGAGACATGGTAAAAAGTCAAAATCTAGATGCAGCAAAAAGCCACTGCATGTCAATTTCAGAGAGCTGGGTTGGGACGATTGGGTCATCGCTCCATTAGATTATGAGGCTTATCACTGTGAGGGCATGTGTGACTTTCCCCTCCGATCTCACCTGGAACCAACCAATCATGCCATCATACAAACTCTAATGAACTCAATGAACCCCAGCAACATGCCACCCAGCTGTTGCGTCCCCTCCAAACTCAGTCCCATTAGCATCTTGTACATTGACGCAGGAAATAATGTTGTGTACAAGCAGTATGAAGACATGGTAGTGGAGTCCTGCGGCTGCAGATGA"

    # Split the sequence into exon1, intron, and exon2
    parts = full_sequence.upper().split('^')
    if len(parts) != 3:
        print("Error: The sequence format with '^' delimiters is incorrect.")
        return

    exon2_seq = parts[2]

    choices = [
        "A. AGCGGTTTACTGAGACCCGG(TGG)", "B. TCCGGCGGGTTTTCGAGTGGG",
        "C. TTCATGCCCCTGGATGCGCT(TGG)", "D. CAGGACCGGTTTCAGATGCG(CGG)",
        "E. GCATCTGAAACCGGTCCTG(TGG)", "F. GGAAGCAATCCTCCGAACGT(TGG)",
        "G. ACGTTGCGAGGACAGAGTCA(AGG)", "H. CCCTTTCACAAATCCTTCCT(TGG)",
        "I. TTCACCCGCACCTTGAACGG(AGG)", "J. CTTTCTTTCTTTCTTTCTTTC(TTT)",
        "K. CTGCTCTACCCAATCGCCA(TGG)", "L. TGCCTG(CGG)",
        "M. TGCAAAGTAGATCGAGATGG(AGG)", "N. ACAGTCCAGAAGGGCGATCA(AGG)",
        "O. ATG(ACC)"
    ]

    print("--- Analyzing potential sgRNA targets in Exon 2 ---")
    valid_candidates = []

    for choice_str in choices:
        match = re.match(r"([A-Z])\.?\s*([ATGC]+)\s*\(([ATGC]+)\)", choice_str)
        if not match:
            # Handle cases without parentheses like choice B
            if choice_str.startswith("B."):
                 match = re.match(r"([A-Z])\.?\s*([ATGC]{17,20})([ATGC]{3})", choice_str)
                 if not match: continue
            else: continue

        info = {
            'id': match.group(1),
            'guide': match.group(2),
            'pam': match.group(3)
        }

        # Rule 1: Must be an NGG PAM
        if not (info['pam'].endswith('GG')):
            continue

        target_seq = info['guide'] + info['pam']
        
        # Rule 2: Must be present in exon 2
        position = exon2_seq.find(target_seq)
        if position != -1:
            # It's a valid candidate. Now analyze its quality.
            guide_seq = info['guide']
            gc_count = guide_seq.count('G') + guide_seq.count('C')
            gc_percent = (gc_count / len(guide_seq)) * 100 if len(guide_seq) > 0 else 0
            
            # Rule 3: Avoid poly-T stretches (which can terminate guide transcription)
            has_poly_t = 'TTTT' in guide_seq

            # Rule 4: GC content should ideally be 40-80%
            good_gc = 40 <= gc_percent <= 80

            valid_candidates.append({
                'choice_str': f"{info['id']}. {info['guide']}({info['pam']})",
                'id': info['id'],
                'position': position,
                'gc_percent': gc_percent,
                'has_poly_t': has_poly_t,
                'good_gc': good_gc,
                'guide': info['guide'],
                'pam': info['pam'],
                'target_seq': target_seq
            })
    
    print("\n--- Evaluation of On-Target Hits Found in Exon 2 ---")
    if not valid_candidates:
        print("No valid targets from the choices were found in Exon 2.")
        return

    # Sort candidates by position (earlier in exon is better)
    valid_candidates.sort(key=lambda x: x['position'])
    
    for candidate in valid_candidates:
        print(f"Candidate {candidate['choice_str']}:")
        print(f"  - Found at position: {candidate['position']}")
        print(f"  - GC Content: {candidate['gc_percent']:.1f}% ({'Good' if candidate['good_gc'] else 'Suboptimal'})")
        print(f"  - Contains 'TTTT' tract: {'Yes' if candidate['has_poly_t'] else 'No'}")
        print("-" * 20)

    # Select the best candidate (the one found earliest in the exon)
    best_candidate = valid_candidates[0]
    
    print("\n--- Recommendation ---")
    print(f"The most suitable target is choice {best_candidate['id']}.")
    print("Reasoning: It is found earliest in the exon, maximizing the chance of a functional knockout.")
    print("It also has ideal biochemical properties (good GC content, no disruptive poly-T tracts).")

    print("\nFinal Equation:")
    final_eq_str = (
        f"Choice {best_candidate['id']}: "
        f"{best_candidate['guide']} + {best_candidate['pam']} = {best_candidate['target_seq']}"
    )
    print(final_eq_str)


if __name__ == "__main__":
    analyze_sgrna_targets()
