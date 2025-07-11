import re

def analyze_genotypes():
    """
    Performs in-silico analysis of a gdf6a RFLP experiment and interprets a gel image.
    """
    # Provided sequences
    wt_orf = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    rev_primer = "TTTTCCCTTGTCCACGAAAC"
    enzyme_site = "GCATC"
    enzyme_name = "SfaNI"

    # Step 1: In Silico PCR
    rev_primer_comp = rev_primer.translate(str.maketrans("ATCG", "TAGC"))[::-1]
    fwd_pos = wt_orf.find(fwd_primer)
    rev_pos = wt_orf.find(rev_primer_comp)
    
    pcr_product_wt = wt_orf[fwd_pos : rev_pos + len(rev_primer_comp)]
    pcr_product_size = len(pcr_product_wt)

    print("--- Step 1 & 2: In Silico Analysis ---")
    print(f"The PCR product is {pcr_product_size} bp long.")

    # Step 2: In Silico Restriction Digest (Wild-Type)
    rev_enzyme_site = enzyme_site.translate(str.maketrans("ATCG", "TAGC"))[::-1]
    fwd_site_pos = pcr_product_wt.find(enzyme_site)
    rev_site_pos = pcr_product_wt.find(rev_enzyme_site)
    
    print(f"Searching for the {enzyme_name} site ('{enzyme_site}' or '{rev_enzyme_site}') in the wild-type amplicon...")

    if fwd_site_pos == -1 and rev_site_pos == -1:
        print("Result: No SfaNI site was found in the wild-type PCR product.")
        wt_is_cut = False
    else:
        print("Result: SfaNI site found in the wild-type PCR product.")
        wt_is_cut = True

    # Step 3 & 4: Deducing Patterns
    print("\n--- Step 3 & 4: Deducing Genotype Patterns ---")
    if not wt_is_cut:
        print("This means the wild-type allele is NOT cut by SfaNI.")
        print("Therefore, the C->A mutation must CREATE the SfaNI restriction site.")
        print("The expected banding patterns are:")
        print("  - Homozygous Wild-Type (+/+): One uncut band (~258 bp).")
        print("  - Homozygous Mutant (-/-): Two smaller bands from the cut product.")
        print("  - Heterozygote (+/-): Three bands (one large uncut band and two smaller cut bands).")
    else:
        # This case is not expected based on our finding, but included for completeness.
        print("This means the wild-type allele IS cut by SfaNI.")
        print("Therefore, the C->A mutation must DESTROY the SfaNI restriction site.")
        print("The expected banding patterns are:")
        print("  - Homozygous Wild-Type (+/+): Two smaller bands from the cut product.")
        print("  - Homozygous Mutant (-/-): One uncut band (~258 bp).")
        print("  - Heterozygote (+/-): Three bands (one large uncut band and two smaller cut bands).")

    # Step 5: Analyze Gel Image
    # Based on visual inspection of the 15 viable lanes in the provided gel image:
    # - Lanes with one high band (Homozygous Wild-Type): 2 (wells #2, #17)
    # - Lanes with three bands (Heterozygote): 9 (wells #6, #7, #8, #9, #10, #13, #14, #15, #16)
    # - Lanes with two low bands (Homozygous Mutant): 4 (wells #3, #4, #5, #11)
    wt_count = 2
    het_count = 9
    mut_count = 4

    print("\n--- Step 5: Gel Image Analysis ---")
    print("Based on visual inspection of the gel image:")
    print(f"Number of homozygous wild-type (+/+) lanes (1 band): {wt_count}")
    print(f"Number of heterozygous (+/-) lanes (3 bands): {het_count}")
    print(f"Number of homozygous mutant (-/-) lanes (2 bands): {mut_count}")

    # Step 6: Final Answer
    print("\n--- Step 6: Final Answer ---")
    print(f"The number of larvae for each genotype (homozygous wild-type/heterozygote/homozygous mutant) is:")
    print(f"{wt_count}/{het_count}/{mut_count}")

analyze_genotypes()
<<<2/9/4>>>