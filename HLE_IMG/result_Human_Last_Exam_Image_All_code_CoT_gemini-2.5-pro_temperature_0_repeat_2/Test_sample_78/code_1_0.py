def analyze_genotypes():
    """
    Analyzes zebrafish gdf6a sequences to predict restriction digest patterns
    and interprets a gel image to determine progeny genotypes.
    """
    # Define the sequences and primers from the problem description.
    wt_orf_provided = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    fwd_primer = "TTTTACGCGCTCTTCGTTTT"
    rev_primer_seq = "TTTTCCCTTGTCCACGAAAC"
    
    # Reverse complement the reverse primer to find it in the ORF.
    rev_primer_rc = rev_primer_seq.translate(str.maketrans("ATCG", "TAGC"))[::-1]

    print("--- Step 1: In Silico PCR and Restriction Analysis ---")
    
    # The problem states the WT gene has a 'C' at position 164, but the provided sequence has a 'G'.
    # This is likely a typo in the provided sequence string. We will proceed assuming the text description is correct.
    mutation_pos_1based = 164
    mutation_pos_0based = mutation_pos_1based - 1
    
    # Create the corrected wild-type sequence and the mutant sequence.
    wt_orf_corrected = wt_orf_provided[:mutation_pos_0based] + 'C' + wt_orf_provided[mutation_pos_0based+1:]
    
    # Find the amplicon.
    fwd_start = wt_orf_corrected.find(fwd_primer)
    rev_start = wt_orf_corrected.find(rev_primer_rc)
    wt_amplicon = wt_orf_corrected[fwd_start : rev_start + len(rev_primer_rc)]
    pcr_product_size = len(wt_amplicon)
    
    print(f"The PCR product is predicted to be {pcr_product_size} bp long.")
    
    # Analyze for the SfaNI site and calculate fragment sizes.
    sfaNI_site = "GCATC"
    wt_site_pos = wt_amplicon.find(sfaNI_site)
    
    # SfaNI cuts as GCATC(N)5/, meaning it cuts 5 bases downstream of the recognition sequence.
    cut_site_offset = len(sfaNI_site) + 5
    cut_position_in_amplicon = wt_site_pos + cut_site_offset
    fragment1_size = cut_position_in_amplicon
    fragment2_size = pcr_product_size - fragment1_size
    
    print(f"The C->A mutation at position 164 destroys the SfaNI recognition site ({sfaNI_site}).")
    print("Therefore, the wild-type allele will be cut, while the mutant allele will not.")
    print(f"Digestion of the wild-type allele produces two fragments: {fragment1_size} bp and {fragment2_size} bp.")

    print("\n--- Step 2: Predicted Gel Patterns ---")
    print(f"1. Homozygous Wild-Type (+/+): Two bands at {fragment1_size} bp and {fragment2_size} bp.")
    print(f"2. Homozygous Mutant (-/-): One uncut band at {pcr_product_size} bp.")
    print(f"3. Heterozygous (+/-): Three bands at {pcr_product_size} bp, {fragment1_size} bp, and {fragment2_size} bp.")

    print("\n--- Step 3: Gel Interpretation and Final Counts ---")
    print("By comparing these predictions to the gel image, we can count the individuals for each genotype:")
    
    # Counts are determined by manually inspecting the gel image provided.
    # A = Homozygous wild-type (two low bands)
    # B = Heterozygous (three bands)
    # C = Homozygous mutant (one high band)
    homozygous_wt_count = 4
    heterozygous_count = 10
    homozygous_mutant_count = 2
    
    print(f"Number of homozygous wild-type larvae (A): {homozygous_wt_count}")
    print(f"Number of heterozygous larvae (B): {heterozygous_count}")
    print(f"Number of homozygous mutant larvae (C): {homozygous_mutant_count}")
    
    print("\nThe final answer in the format A/B/C is:")
    print(f"{homozygous_wt_count}/{heterozygous_count}/{homozygous_mutant_count}")

analyze_genotypes()
<<<4/10/2>>>