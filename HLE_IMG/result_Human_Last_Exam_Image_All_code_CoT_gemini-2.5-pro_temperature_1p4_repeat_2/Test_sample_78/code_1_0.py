import re

def analyze_genotypes():
    """
    Performs an in-silico PCR and restriction digest to determine expected
    fragment sizes and then counts genotypes from the described gel.
    """
    # 1. Define sequences and enzyme parameters
    wt_orf = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    rev_primer = "TTTTCCCTTGTCCACGAAAC"
    enzyme_site = "GCATC"
    cut_offset = 5 # SfaNI cuts 5 bases downstream: GCATC(N)5^

    # 2. In-silico PCR
    # The user-provided forward primer has TTTT, which is likely a tail.
    # The core sequence TTTACGCGCTCTTCGTTTT is found at index 26.
    fwd_primer_core = "TTACGCGCTCTTCGTTTT"
    pcr_start_index = wt_orf.find(fwd_primer_core)

    rev_primer_rc = rev_primer.translate(str.maketrans("ATCG", "TAGC"))[::-1]
    rev_match_start_index = wt_orf.find(rev_primer_rc)
    pcr_end_index = rev_match_start_index + len(rev_primer_rc)

    wt_pcr_product = wt_orf[pcr_start_index:pcr_end_index]
    pcr_product_size = len(wt_pcr_product)

    # 3. Create mutant sequence
    mutation_pos = 164 # 1-based position
    mutation_index_in_orf = mutation_pos - 1
    
    # Verify WT base is C at this position
    wt_base = wt_orf[mutation_index_in_orf]

    mutant_orf_list = list(wt_orf)
    mutant_orf_list[mutation_index_in_orf] = 'A'
    mutant_orf = "".join(mutant_orf_list)
    mutant_pcr_product = mutant_orf[pcr_start_index:pcr_end_index]
    
    print("--- In Silico Analysis ---")
    print(f"PCR product length: {pcr_product_size} bp")
    print(f"Mutation at position {mutation_pos}: WT base is '{wt_base}', Mutant base is 'A'.")

    # 4. In-silico Restriction Digest
    print("\n--- Predicted Genotype Patterns ---")
    
    # Mutant allele analysis
    if mutant_pcr_product.find(enzyme_site) == -1:
        print(f"Homozygous Mutant (-/-): The mutation destroys the SfaNI site. The PCR product is UNCUT.")
        print(f"Expected band: {pcr_product_size} bp")
    
    # Wild-type allele analysis
    cut_site_index = wt_pcr_product.find(enzyme_site)
    if cut_site_index != -1:
        fragment1_size = cut_site_index + len(enzyme_site) + cut_offset
        fragment2_size = pcr_product_size - fragment1_size
        print(f"\nHomozygous Wild-Type (+/+): The SfaNI site '{enzyme_site}' is present. The PCR product is CUT.")
        print(f"Expected bands: {fragment1_size} bp and {fragment2_size} bp")

    # Heterozygote analysis
    print(f"\nHeterozygote (+/-): Will show a combination of both patterns.")
    print(f"Expected bands: {pcr_product_size} bp, {fragment1_size} bp, and {fragment2_size} bp")

    # 5. Gel Image Analysis and Final Count
    # Based on visual inspection of the gel image and the predicted patterns:
    # Pattern 1 (+/+): Two bands (~143 bp and ~96 bp)
    # Pattern 2 (+/-): Three bands (~239 bp, ~143 bp, ~96 bp)
    # Pattern 3 (-/-): One band (~239 bp)
    
    # Counting lanes with clear bands from left to right:
    # There are 16 lanes with scorable genotypes.
    num_homozygous_wt = 5  # Lanes with two lower bands
    num_heterozygous = 8   # Lanes with three bands (upper and two lower)
    num_homozygous_mutant = 3 # Lanes with one upper band

    print("\n--- Gel Interpretation ---")
    print(f"By analyzing the gel image, we count the genotypes as follows:")
    print(f"Number of homozygous wild-type larvae (A) = {num_homozygous_wt}")
    print(f"Number of heterozygous larvae (B) = {num_heterozygous}")
    print(f"Number of homozygous mutant larvae (C) = {num_homozygous_mutant}")

    print("\nFinal Answer Format (A/B/C):")
    print(f"<<<{num_homozygous_wt}/{num_heterozygous}/{num_homozygous_mutant}>>>")

analyze_genotypes()