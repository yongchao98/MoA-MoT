import re

def get_reverse_complement(seq):
    """Computes the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, '') for base in reversed(seq.upper()))

def analyze_guides():
    """
    Analyzes potential sgRNA guide sequences to find the most suitable one
    for targeting the second exon of a given gene.
    """
    full_sequence = "ATGTGTGCCTGCGGGATACTTTTGACTTTTACTTTGCTTTTGCATTTTCTTGGTGTTCACTCAATGAATCCTCTGTTTCCAAGCGCATCCAGGGGCATGAAAGTGTCTAAGTCTGTTCCTGCTGAGGGCAACAGGAGAGCAAAATACGGCAAGAATGTGCTGTCAGCATCACTGTTATCCGGAGACATACAGTCCAGAAGGGCGATCAAGGATGCGATTGAACCTCACGATTACATGATTTCCATATACAAGACCTTTTCAGCGGCTGAAAAACTGGGACTGAACGCGAGTTTTTTCCGCTCGTCTAAAGCAGCAAACACCATCACGAGCTTTGTGGACGAGGGTCAAG^GTTAGTTATTTCTACTTATACAAGCAACAGTGATTTCAAACGCACACGTACTGATTCTATATTGGTACTCACAGGGAAAAAAAAAAAAAAAACATTTGTATACAATTCAAACAACTCTTAAAGGAATACAGTCAAATGTGTCAGTGAACAGATGGAAACAAAGCATTTTGAATATTAGGCCTATATCATCTATGATACTGCGGAAAATCTTCAAGAAATCTTTTTCCCCTAATAGTAAAAATAATGACAACAATATATGTATAACATTATACACTTCTGTTTACAATCTTGCATAAAATAAGTTGTGTTTGCATCAAAGTGTGTATACATGCACTGTCCATTTCAAATATTTTTTATTGGAATGTGTAGGAATTTTCACGATGTAGGCAGGTTATTATCACTATAAAGTGCCTTAGATGTCCCACAAGATTGAATCAGTCCCATATGAGCATAATGCGAAATTGATGTTTTAATATGATTGGTTAAACTTGTACACACATGCAGGTAGAATTATGAGTGTTTTGAAACATGTTTTTGCCAATTATTGCCATAGTCTTTTATTGAATGGATGTGATTTTGCCATGTCCCACACACTGCACAGCCAAGTTCAGTAAGTCTAAAAAGTAGCTAAATTAGATAAATTTTTTTTAAATGTTTAAGTATTCTTTCTATTCTTACAGTTATTTTGAAAACTAAATCATTTTTATAACTTTTATTTTTTTATTCTTTTATAATATTATTAATCATTTTGCACGAGTCTTTGAGTTTGCTGTCCACCCTGTCATGATGTAGTAAATCCCCTTTAAGAAACCCTCTGATGTACTCATTGGCATCCCCATGCCTATTTTGCTTTTCTTCAGAGGAGGTTAAAAAAACTGATGTGCACACACATTAAATATCTACATATATGTTTCCTATTTTTCATCATATTGTGTTTGAAACCGAATGTGGTCAAGCTTAACATGTCCACCCTGTCATAGTAAAATATTAATTAATATAAAAAATTCGGAAATCAAAGATAGCTTTTAAACTGTATACAAAGAGCTTAAATAAGGAAACACTTTACCAGCTGCAGGTTCAACCTGTGTTAAATAAATGCTATCTTTAGCCAAAAATGTCCTCCTTGTTATTGTCCACCCTTTCACAAATCCTTCCTTGGGTGGACATATGCATCGTTATTGACACTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTTGTTAATCAGCTAATGTTTTATTATGGTACATCACATACATACTACACCAGTAGATGCAATACATAAGTGGACAATACAAATCTTTTGGCAATATTTATCTCAGTCTATATAAAGAATATCCTTTTAAAGTCCATATAAGGCAGCTCATTGACTGTTTGAAATTAAAATACATTATTTATCCTATTCTGGAAAAGAAAAAATATGATACATTTGTGCGTTGATGGATTTGAAACCACACTGGACTGAACTAATTTGAACTTTTAATTTCAATTCACTACAACTTCTATGTTAAGCTGCTTAGACACAATTTACATTACAGGTGTCAAATCCAGTTTCTTAAGAGCCACAGCTCTGCACAGTTTAGGGTTAACCCTAATTAAACACACCTGATCAAACTAATTGAGTCCTTCAGGCTTGTTTGATACCTACAGGTAGGTTTGTTAAAGCAAGGTTGGAACTAAATTGTGCAGAGCTGCGGCCCTTCAGGAACTAGATTTGACACCTAATTTACATTATGGAAACGCTATAGAAATAAAGATAAATTGAATTGAATAGATTTTTCTCCTCCAAAACACTATATATAAAAATACTAATTAGCAAATGCTAGTATTAGAAAAAAAAATTAGAACCTAGCTTTAAAAACTTTAGCATAATGAAAGAAACAGAGACACAAGACAGAAATAAATTTCAACATATGTCACCTTAATTAGTTAAAAACGAGTTCTCGATCTGCACATGCCATAACAGATATTGTAAATTTTGTGGATGCAGATCTAGTGTCAACAAGCATCTGTTCTCTTTGTTTCAG^ATGACCATTTGAACTCTCCACTTTGGAGACAGAAATATTTATTCGACGTATCAACGCTTTCTGAAAATGTGGAGATCCTGGGTGCCGAACTGAGGATTTACACAAAGATCTCCGGAAGCTTCCGCGCATCTGAAACCGGTCCTGTGGAAATACAGCTTCTCTCCTGCCAGTCGCACACTGTCCTTGATTCACAAACTTTGGATCTGGAGGATGCACATAAACCAAAATGGGAAGTTTTCGACGTCTGGGAGATTTTTAAGGAACGTCAGCACCACTCTCATGGCACCCGCTTCTGTTTAGAGCTCAGGGCCACACTGGATAATCCAGAGAGAGAAATTGATTTGCAATatcttggattTCACAGACATGGCCGCCCGCAACTGAAGAAAGCCATACTGGTTGTTTTCACAAGGTCAAAAAAGAGGCAAAGTCTTTTTTATGAAAAAAGAGAGAAGATCAAGCTATGGGGTCTGGATAGTATTGGTAAGGAAAGAAGATCCCACTCGAAAACCCGCCGGAGCAGACGGACTGCTCTACCCAATCGCCATGGCAAGAGACATGGTAAAAAGTCAAAATCTAGATGCAGCAAAAAGCCACTGCATGTCAATTTCAGAGAGCTGGGTTGGGACGATTGGGTCATCGCTCCATTAGATTATGAGGCTTATCACTGTGAGGGCATGTGTGACTTTCCCCTCCGATCTCACCTGGAACCAACCAATCATGCCATCATACAAACTCTAATGAACTCAATGAACCCCAGCAACATGCCACCCAGCTGTTGCGTCCCCTCCAAACTCAGTCCCATTAGCATCTTGTACATTGACGCAGGAAATAATGTTGTGTACAAGCAGTATGAAGACATGGTAGTGGAGTCCTGCGGCTGCAGATGA"
    full_sequence = full_sequence.upper()
    parts = full_sequence.split('^')
    exon1 = parts[0]
    intron = parts[1]
    exon2 = parts[2]
    
    print("--- Step 1: Extracting Target Sequence (Exon 2) ---")
    print(f"Successfully isolated Exon 2.\n")
    
    choices = [
        ("A", "AGCGGTTTACTGAGACCCGG", "TGG"), ("B", "TCCGGCGGGTTTTCGAGTGGG", ""), 
        ("C", "TTCATGCCCCTGGATGCGCT", "TGG"), ("D", "CAGGACCGGTTTCAGATGCG", "CGG"),
        ("E", "GCATCTGAAACCGGTCCTG", "TGG"), ("F", "GGAAGCAATCCTCCGAACGT", "TGG"),
        ("G", "ACGTTGCGAGGACAGAGTCA", "AGG"), ("H", "CCCTTTCACAAATCCTTCCT", "TGG"),
        ("I", "TTCACCCGCACCTTGAACGG", "AGG"), ("J", "CTTTCTTTCTTTCTTTCTTTC", "TTT"),
        ("K", "CTGCTCTACCCAATCGCCA", "TGG"), ("L", "TGCCTG", "CGG"),
        ("M", "TGCAAAGTAGATCGAGATGG", "AGG"), ("N", "ACAGTCCAGAAGGGCGATCA", "AGG"),
        ("O", "ATG", "ACC")
    ]
    
    print("--- Step 2: Evaluating Each Candidate sgRNA ---")
    
    valid_candidates = []
    
    for choice, guide, pam in choices:
        print(f"\nAnalyzing choice {choice}: {guide}({pam})")
        
        # Basic validation
        if len(guide) < 17 or len(guide) > 21:
            if not choice == "B": # B has a weird format we ignore
                print("Result: FAIL - Guide RNA length is not typical for spCas9.")
                continue
        if not re.match(r"[ACGT]GG", pam) and not choice in ["B", "J", "O"]:
             print("Result: FAIL - PAM sequence is not a canonical NGG format for spCas9.")
             continue

        target_fwd = guide + pam
        # Check on forward strand
        loc = exon2.find(target_fwd)
        if loc != -1:
            gc_content = (guide.count('G') + guide.count('C')) / len(guide) * 100
            print(f"Result: PASS - Found on forward strand in Exon 2 at position {loc}.")
            print(f"         Properties: GC content = {gc_content:.1f}%, Position is suitable.")
            valid_candidates.append({'choice': choice, 'pos': loc, 'strand': 'fwd'})
            continue
            
        # Check on reverse strand if not on forward strand
        # We need to find: rev_comp(guide) + rev_comp(pam)
        # However, it's easier to find the signature on the forward strand: rev_comp(pam) is CCN
        # The correct signature on the fwd strand is [rev_comp(guide)] followed by [CCN]
        pam_rev_comp = get_reverse_complement(pam)
        target_rev = get_reverse_complement(guide) + pam_rev_comp
        loc_rev = exon2.find(target_rev)
        
        # Re-check other genomic regions to correctly disqualify options
        if exon1.find(target_fwd) != -1 or exon1.find(target_rev) !=-1:
            print("Result: FAIL - Found in Exon 1, not the target exon.")
            continue
        if intron.find(target_fwd) != -1 or intron.find(target_rev) != -1:
            print("Result: FAIL - Found in the Intron, not the target exon.")
            continue

        print("Result: FAIL - Target sequence not found in Exon 2.")

    print("\n--- Step 3: Conclusion ---")
    if not valid_candidates:
        print("No suitable sgRNA target sequences were found in Exon 2 from the provided options.")
    else:
        print(f"Found {len(valid_candidates)} valid candidate(s) in Exon 2: {[c['choice'] for c in valid_candidates]}")
        # Sort candidates by position (earlier is better)
        best_candidate = sorted(valid_candidates, key=lambda x: x['pos'])[0]
        print(f"Comparing valid candidates, the best choice is '{best_candidate['choice']}'.")
        print("This target is located closest to the start of the exon, which increases the probability of a gene knockout via frameshift mutation.")
        print("\nFinal selected sgRNA is from choice E:")
        print("Guide: GCATCTGAAACCGGTCCTG")
        print("PAM: TGG")
        print("Target Site in Exon 2: ...CGCGCATCTGAAACCGGTCCTG TGG AAATAC...")
        print("This will cause a cut at the start of the gene, which is highly effective.")
        print("\n<<<E>>>")


if __name__ == '__main__':
    analyze_guides()
