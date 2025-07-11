def analyze_genotypes():
    """
    Performs an in silico PCR and restriction digest to determine the expected
    banding patterns and then reports the count of genotypes from the gel image.
    """
    # 1. In Silico Analysis
    
    # Provided sequences and information
    orf_wt = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    rev_primer = "TTTTCCCTTGTCCACGAAAC"
    
    # Find PCR product
    fwd_pos = orf_wt.find(fwd_primer)
    rev_primer_rc = "GTTTCGTGGACAAGGGAAAA" # Reverse complement of rev_primer
    rev_pos = orf_wt.find(rev_primer_rc)
    
    pcr_product_wt = orf_wt[fwd_pos : rev_pos + len(rev_primer_rc)]
    pcr_product_size = len(pcr_product_wt)
    
    # Create mutant sequence
    mutation_pos_orf = 164  # 1-based position in ORF
    mutation_pos_pcr = mutation_pos_orf - (fwd_pos + 1) # 0-based index in PCR product
    
    if pcr_product_wt[mutation_pos_pcr] == 'C':
        pcr_product_mut = pcr_product_wt[:mutation_pos_pcr] + 'A' + pcr_product_wt[mutation_pos_pcr+1:]
    else:
        print("Error: Expected 'C' at mutation site but found another base.")
        return

    # In silico restriction digest with SfaNI (Site: GCATC, Cleavage: GCATC(N)5^)
    sfaNI_site = "GCATC"
    site_pos_wt = pcr_product_wt.find(sfaNI_site)
    
    print("--- In Silico Analysis ---")
    
    # Wild-type allele
    if site_pos_wt != -1:
        # Cleavage is after N5, so 5 bases after the 5-base recognition site
        cut_site = site_pos_wt + 5 + 5
        fragment1_size = cut_site
        fragment2_size = pcr_product_size - fragment1_size
        print(f"Wild-Type Allele: Contains an SfaNI site. The {pcr_product_size} bp product is cut into {fragment1_size} bp and {fragment2_size} bp fragments.")
    else:
        print(f"Wild-Type Allele: No SfaNI site found. Product is uncut at {pcr_product_size} bp.")
        fragment1_size, fragment2_size = None, None

    # Mutant allele
    site_pos_mut = pcr_product_mut.find(sfaNI_site)
    if site_pos_mut == -1:
        print(f"Mutant Allele: The C->A mutation destroys the SfaNI site. Product is uncut at {pcr_product_size} bp.")
    else:
        print("Mutant Allele: SfaNI site is unexpectedly present.")

    # 2. Predicted Genotype Patterns
    print("\n--- Predicted Gel Patterns ---")
    print(f"Homozygous Wild-Type (WT/WT): Two bands at {fragment1_size} bp and {fragment2_size} bp.")
    print(f"Homozygous Mutant (MUT/MUT): One band at {pcr_product_size} bp (uncut).")
    print(f"Heterozygote (WT/MUT): Three bands at {pcr_product_size} bp, {fragment1_size} bp, and {fragment2_size} bp.")
    
    # 3. Gel Interpretation and Count
    # Based on careful observation of the gel image provided.
    # Pattern 1 (WT/WT): Two low bands (~147 bp, ~79 bp)
    # Pattern 2 (WT/MUT): Three bands (~226 bp, ~147 bp, ~79 bp)
    # Pattern 3 (MUT/MUT): One high band (~226 bp)
    
    # Lane-by-lane analysis:
    # 1: MUT/MUT
    # 2: WT/WT
    # 3: WT/MUT
    # 4: WT/WT
    # 5: WT/WT
    # 6: WT/MUT
    # 7: WT/MUT
    # 8: WT/MUT
    # 9: WT/MUT
    # 10: WT/WT
    # 11: WT/MUT
    # 12: WT/MUT
    # 13: WT/MUT
    # 14: WT/WT
    # 15: WT/MUT
    
    num_homozygous_wt = 5
    num_heterozygous = 9
    num_homozygous_mutant = 1
    
    print("\n--- Genotype Count from Gel ---")
    print(f"Number of homozygous wild-type larvae (WT/WT): {num_homozygous_wt}")
    print(f"Number of heterozygous larvae (WT/MUT): {num_heterozygous}")
    print(f"Number of homozygous mutant larvae (MUT/MUT): {num_homozygous_mutant}")
    
    print("\nFinal Answer Format (WT/Het/Mut):")
    print(f"{num_homozygous_wt}/{num_heterozygous}/{num_homozygous_mutant}")

if __name__ == '__main__':
    analyze_genotypes()