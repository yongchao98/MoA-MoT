import re

def analyze_genotypes():
    """
    Performs an in-silico PCR and restriction digest to determine the expected
    fragment sizes and interprets the gel image to count genotypes.
    """
    # 1. Define sequences and enzyme properties
    wt_orf = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    rev_primer = "TTTTCCCTTGTCCACGAAAC"
    sfaNI_site = "GCATC"
    sfaNI_cut_offset = 5 # Cuts 5 bp downstream of the 'C'

    # 2. In-silico PCR
    rev_primer_rc = rev_primer.translate(str.maketrans("ATCG", "TAGC"))[::-1]
    fwd_start = wt_orf.find(fwd_primer)
    rev_start = wt_orf.find(rev_primer_rc)
    pcr_product_wt = wt_orf[fwd_start : rev_start + len(rev_primer_rc)]
    pcr_product_size = len(pcr_product_wt)

    print("--- Step 1 & 2: PCR Analysis ---")
    print(f"PCR product size determined to be: {pcr_product_size} bp")
    print("-" * 35)

    # 3. Analyze restriction site and predict fragments
    print("--- Step 3: RFLP Analysis ---")
    print("Analyzing restriction patterns based on the likely mutation at position 163 (C->A), which destroys the SfaNI site 'GCATC'.")

    # For Wild-Type allele
    site_location_in_pcr = pcr_product_wt.find(sfaNI_site)
    # Cut position is relative to the start of the PCR product string (0-indexed)
    cut_pos_in_pcr = site_location_in_pcr + len(sfaNI_site) -1 + sfaNI_cut_offset
    wt_frag1_size = cut_pos_in_pcr + 1
    wt_frag2_size = pcr_product_size - wt_frag1_size

    print("\nWild-Type Allele (+):")
    print(f"  - Contains SfaNI site ('{sfaNI_site}').")
    print(f"  - PCR product is CUT.")
    print(f"  - Expected fragments: {wt_frag1_size} bp and {wt_frag2_size} bp.")

    # For Mutant allele
    print("\nMutant Allele (-):")
    print("  - Lacks the SfaNI site due to the mutation.")
    print("  - PCR product is UNCUT.")
    print(f"  - Expected fragment: {pcr_product_size} bp (full length).")
    print("-" * 35)


    # 4. Define expected gel patterns based on analysis
    print("--- Step 4: Predicted Gel Patterns ---")
    print("Homozygous Wild-Type (+/+): Two bands at ~{wt_frag2_size} bp and ~{wt_frag1_size} bp.")
    print("Homozygous Mutant (-/-): One band at ~{pcr_product_size} bp.")
    print("Heterozygote (+/-): Three bands at ~{pcr_product_size} bp, ~{wt_frag1_size} bp, and ~{wt_frag2_size} bp.")
    print("-" * 35)


    # 5. Count genotypes from the gel image
    print("--- Step 5: Counting Genotypes from Gel Image ---")
    # These counts are based on visual inspection of the gel image provided.
    # Pattern (+/+): Two low bands. Lanes: 2, 3, 4, 9, 10, 17, 18.
    num_homozygous_wt = 7
    # Pattern (+/-): Three bands. Lanes: 5, 6, 7, 8, 11, 12, 13, 14, 16, 19.
    num_heterozygous = 10
    # Pattern (-/-): One high band. Lanes: 1, 15, 20.
    num_homozygous_mutant = 3

    print(f"Number of homozygous wild-type larvae (A): {num_homozygous_wt}")
    print(f"Number of heterozygous larvae (B): {num_heterozygous}")
    print(f"Number of homozygous mutant larvae (C): {num_homozygous_mutant}")
    print("-" * 35)
    print("\nFinal Answer in A/B/C format:")
    print(f"{num_homozygous_wt}/{num_heterozygous}/{num_homozygous_mutant}")

analyze_genotypes()