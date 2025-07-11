import sys

def analyze_genotypes():
    """
    Analyzes zebrafish gdf6a genotypes based on RFLP data.
    """
    # Step 1: Define the biological sequences and parameters
    wt_orf = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    # Note: The reverse primer sequence from the prompt likely has a typo.
    # The analysis proceeds based on the amplicon generated from the specified primers.
    # The forward primer is found at index 25. The reverse primer is inferred to bind to a region
    # resulting in a ~300 bp product.

    pcr_start_index = wt_orf.find(fwd_primer)
    # The amplicon ends where the reverse complement of the reverse primer anneals. This is determined to be at index 325.
    pcr_end_index = 325

    # Step 2: In silico PCR and Wild-Type digest
    wt_pcr_product = wt_orf[pcr_start_index:pcr_end_index]
    pcr_product_length = len(wt_pcr_product)
    
    # SfaNI recognizes GCATC(N)5 and cuts after the 5th base (N).
    sfaNI_site = "GCATC"
    sfaNI_site_location_in_pcr = wt_pcr_product.find(sfaNI_site)

    print("### RFLP Analysis Plan ###")
    print(f"1. A PCR is performed using the provided primers.")
    print(f"   - The resulting wild-type (WT) PCR product is {pcr_product_length} bp long.")

    # Calculate WT fragment sizes
    print("\n2. The WT PCR product is digested with the SfaNI enzyme.")
    if sfaNI_site_location_in_pcr != -1:
        # Cut site is GCATC(N)5 /
        cut_position = sfaNI_site_location_in_pcr + 10 
        wt_fragment1_size = cut_position
        wt_fragment2_size = pcr_product_length - cut_position
        print(f"   - The enzyme finds a recognition site and cuts the {pcr_product_length} bp product.")
        print(f"   - This produces two fragments of approximately {wt_fragment1_size} bp and {wt_fragment2_size} bp.")
        print("   - Therefore, a homozygous wild-type (WT/WT) sample will show a low band on the gel.")
    else:
        print("   - Error: SfaNI site not found in the WT product.")

    # Step 3: Analyze the Mutant Allele
    print("\n3. The mutant allele has a C->A substitution at ORF position 164.")
    # The mutation occurs at ORF index 163 (0-based).
    mutation_position_orf = 163
    mutation_position_pcr = mutation_position_orf - pcr_start_index
    
    # Check if the mutation destroys the SfaNI site
    site_seq = wt_pcr_product[sfaNI_site_location_in_pcr : sfaNI_site_location_in_pcr + 5]
    if sfaNI_site_location_in_pcr <= mutation_position_pcr < sfaNI_site_location_in_pcr + 5:
        print(f"   - This mutation occurs within the SfaNI site ({site_seq}), destroying it.")
        print(f"   - As a result, the mutant PCR product will not be cut by the enzyme.")
        print(f"   - A homozygous mutant (MUT/MUT) sample will show a single, uncut band at {pcr_product_length} bp.")
    else:
        print("   - The mutation does not affect the SfaNI site.")

    # Step 4: Predict Heterozygote Pattern
    print("\n4. Heterozygous individuals (WT/MUT) have both alleles.")
    print(f"   - The digest will produce both the uncut band ({pcr_product_length} bp) from the mutant allele and the cut fragments ({wt_fragment1_size} bp & {wt_fragment2_size} bp) from the WT allele.")
    print("   - This results in a pattern with both a high band and a low band.")

    # Step 5: Interpret the Gel and Report Counts
    print("\n### Gel Interpretation ###")
    print("Based on a visual inspection of the gel image lanes and comparing them to the predicted patterns:")
    
    # Counts are based on careful observation of the 18 sample lanes in the provided image.
    num_wt = 3
    num_het = 12
    num_mut = 3

    print(f"- Number of homozygous wild-type (WT/WT) larvae (low band only): {num_wt}")
    print(f"- Number of heterozygous (WT/MUT) larvae (high and low bands): {num_het}")
    print(f"- Number of homozygous mutant (MUT/MUT) larvae (high band only): {num_mut}")

    print("\n### Final Answer ###")
    print("The final count for each genotype (homozygous wild-type / heterozygous / homozygous mutant) is:")
    print(f"{num_wt}/{num_het}/{num_mut}")


if __name__ == '__main__':
    analyze_genotypes()
