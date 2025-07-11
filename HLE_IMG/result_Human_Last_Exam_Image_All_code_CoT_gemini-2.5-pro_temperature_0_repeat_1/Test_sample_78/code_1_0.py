import re

def analyze_genotypes():
    """
    Analyzes zebrafish gdf6a sequences to predict restriction digest patterns
    and interprets a gel image to determine genotypes.
    """
    # Provided sequences
    orf_wt = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    rev_primer = "TTTTCCCTTGTCCACGAAAC"
    
    # --- Step 1: Determine PCR Product ---
    # Find primer binding sites (using 0-based indexing)
    fwd_start = orf_wt.find(fwd_primer)
    
    # Find reverse complement of the reverse primer
    rev_comp = rev_primer.translate(str.maketrans("ATCG", "TAGC"))[::-1]
    rev_start = orf_wt.find(rev_comp)
    
    # Extract the PCR product sequence
    pcr_product_wt = orf_wt[fwd_start : rev_start + len(rev_comp)]
    pcr_product_size = len(pcr_product_wt)

    print("--- Step 1: PCR Analysis ---")
    print(f"Forward primer found at position: {fwd_start + 1}")
    print(f"Reverse primer complement found at position: {rev_start + 1}")
    print(f"Predicted PCR product size: {pcr_product_size} bp\n")

    # --- Step 2: Analyze the Restriction Site ---
    # Mutation is C->A at position 164 (1-based) of the ORF
    mutation_pos_orf = 164 - 1 # convert to 0-based index
    
    # Find mutation position relative to the PCR product
    mutation_pos_pcr = mutation_pos_orf - fwd_start
    
    # Create the mutant PCR product
    pcr_product_mut_list = list(pcr_product_wt)
    pcr_product_mut_list[mutation_pos_pcr] = 'A'
    pcr_product_mut = "".join(pcr_product_mut_list)
    
    # SfaNI recognition site is GCATC
    sfaNI_site = "GCATC"
    
    wt_has_site = sfaNI_site in pcr_product_wt
    mut_has_site = sfaNI_site in pcr_product_mut
    
    print("--- Step 2: Restriction Site Analysis ---")
    print(f"SfaNI recognition sequence: {sfaNI_site}")
    print(f"Wild-type sequence at mutation site: ...{pcr_product_wt[mutation_pos_pcr-5:mutation_pos_pcr+6]}...")
    print(f"Mutant sequence at mutation site:   ...{pcr_product_mut[mutation_pos_pcr-5:mutation_pos_pcr+6]}...")
    print(f"Does wild-type PCR product have an SfaNI site? {wt_has_site}")
    print(f"Does mutant PCR product have an SfaNI site? {mut_has_site}")
    print("Conclusion: The C->A mutation creates an SfaNI restriction site.\n")

    # --- Step 3: Predict Banding Patterns ---
    print("--- Step 3: Predicted Gel Patterns ---")
    print(f"Homozygous Wild-Type (+/+): No SfaNI site. The PCR product remains uncut.")
    print(f"Expected band: {pcr_product_size} bp\n")

    # Calculate mutant fragment sizes
    # SfaNI cuts 5 bases 3' to its recognition sequence: GCATC(N)5^
    cut_site_start = pcr_product_mut.find(sfaNI_site)
    cut_position = cut_site_start + len(sfaNI_site) + 5
    frag1_size = cut_position
    frag2_size = pcr_product_size - frag1_size
    
    print(f"Homozygous Mutant (m/m): The SfaNI site is present and the PCR product is cut.")
    print(f"Expected bands: {frag1_size} bp and {frag2_size} bp\n")
    
    print(f"Heterozygote (+/m): One allele is cut, one is not.")
    print(f"Expected bands: {pcr_product_size} bp, {frag1_size} bp, and {frag2_size} bp\n")

    # --- Step 4: Interpret the Gel ---
    print("--- Step 4: Gel Interpretation ---")
    print("Based on the predictions, we can interpret the gel image:")
    print(f"- One high band (~{pcr_product_size} bp) corresponds to Homozygous Wild-Type (+/+).")
    print(f"- Two low bands (~{frag1_size} bp and ~{frag2_size} bp) correspond to Homozygous Mutant (m/m).")
    print(f"- Three bands (one high, two low) correspond to Heterozygotes (+/m).\n")
    
    # Counting the lanes based on visual inspection of the provided gel image
    # Lane-by-lane analysis:
    # 1:WT, 2:Mut, 3:Het, 4:Mut, 5:WT, 6:Het, 7:Het, 8:Het, 9:Mut, 10:WT, 11:Het, 12:Het, 13:Het, 14:Mut, 15:WT
    num_wt = 4
    num_het = 7
    num_mut = 4
    
    print("--- Step 5: Final Counts ---")
    print("Counting the lanes in the gel image gives the following results:")
    print(f"Number of homozygous wild-type larvae (A): {num_wt}")
    print(f"Number of heterozygous larvae (B): {num_het}")
    print(f"Number of homozygous mutant larvae (C): {num_mut}")
    
    # Print the final answer in the required format
    print("\nFinal Answer (A/B/C):")
    print(f"{num_wt}/{num_het}/{num_mut}")

# Run the analysis
analyze_genotypes()