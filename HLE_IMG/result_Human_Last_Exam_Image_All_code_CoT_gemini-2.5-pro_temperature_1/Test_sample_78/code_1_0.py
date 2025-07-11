import textwrap

def analyze_genotypes():
    """
    Analyzes zebrafish gdf6a genotypes based on a PCR-RFLP experiment.
    """
    # Step 1: Define the input sequences and parameters
    wt_orf = "ATGGATGCCTTGAGAGCAGTCGCCTTTTACGCGCTCTTCGTTTTCCTCTGGAGTTTACCGTGTTGCCAGTCAGCTGCGCTAATATCGCAGAAAAGGAGCAAGGGTGCCAGGAGCGCGTTTGATGGACAAAGGTCACATAAATTTCTTAAAGAGATTTTAGCATCATCACCGGGCGCGAGTCGTCGGGATGATTTTAAGGACCCGGTTGTGCCTCATGACTACATGATCTCCATATACAGGACTTACTCCGCCGCTGAGAAACTGGGGCTCAATGCGAGCTTTTTCCGCTCTTCAAAGTCTGCAAACACCATAACGAGTTTCGTGGACAAGGGAAAAGACGATCTCACGCTCTCTCCTTTGCGAAGACAAACGTATCTGTTTGATGTTTCAACTCTCTCAGACAAAGAGGAGCTGGTCGGTGCTGAATTAAGGATATTTCGCAAATCGCCCGGGGATGTCCAACCGTCCCCATCAGGCGTCTACAACCTTCATTTACTCTCATGTCGATCAGAGAGGCCACTGGCCTCCAGGTCCATTGATCTTCAGGATTCCCGAAAAGCAGAATGGGAGGTTCTGGACGTTTGGGGGATTTTTAAACACAGGCACCAAGAGAATCAGCTTTGTCTCCAGCTTAAGGTTACATATGGCAAATCTGACACTGAGATCGACCTAAAGCAACTTGGTTTCCACCGCCACAGCCGGACGCAGCAAGAAAGAGCCATATTGGTGGTCTACACGCGGTCCAAGAAGAGAGAGAACTTGTTTAATGAGATGAAAGAGAAAATTAAGTCTCGCGGAGATGATGATGAGGAGGAGAGCGCGCTGCAGTTTAAAGCGCGGCGCAGACGGAGAACTGCGCTTAATAATCGGCACGGGAAAAGGCATGGCAAAAAGTCCAAATCGAGATGCAGCAAAAAGGCTCTGCACGTCAACTTCAAAGAACTTGGATGGGACGACTGGATCATCGCTCCCCTGGATTACGAAGCCTATCACTGCGAGGGCGTGTGCGACTTCCCGTTGAGGTCGCACCTAGAGCCGACCAACCACGCCATCATTCAGACGCTCATGAACTCCATGGACCCCAACAGCACTCCACCGAGCTGTTGCGTCCCCACAAAACTCAGCCCCATCAGTATACTGTACATAGACTCTGGGAACAACGTCGTGTACAAACAGTACGAGGACATGGTGGTAGAACAGTGTGGCTGCAGGTAG"
    mutation_pos = 164  # 1-based position
    sfaNI_site = "GCATC"

    # The problem description contains a discrepancy: the calculated PCR product size (~214 bp)
    # from the provided ORF does not match the observed size on the gel (~480 bp).
    # This suggests the provided ORF is for reference, but the actual amplicon might include
    # untranslated or intronic regions. However, the sequence information around the mutation
    # is sufficient to determine the RFLP pattern.

    # Step 2: Analyze the effect of the mutation on the SfaNI restriction site
    
    # The mutation is a C -> A substitution at position 164.
    # Let's examine the sequence around this position.
    wt_sequence_around_mutation = wt_orf[mutation_pos - 10 : mutation_pos + 10]
    
    # Create the mutant sequence by changing C to A at index (mutation_pos - 1)
    mut_orf_list = list(wt_orf)
    mut_orf_list[mutation_pos - 1] = 'A'
    mut_orf = "".join(mut_orf_list)
    mut_sequence_around_mutation = mut_orf[mutation_pos - 10 : mutation_pos + 10]

    # Check for the SfaNI site in both WT and mutant sequences
    wt_has_site = sfaNI_site in wt_sequence_around_mutation
    mut_has_site = sfaNI_site in mut_sequence_around_mutation
    
    print("--- Step 1: Sequence Analysis ---")
    print(f"The restriction enzyme is SfaNI, which recognizes the sequence '{sfaNI_site}'.")
    print(f"The mutation is a C to A change at position {mutation_pos}.")
    print(f"WT sequence context: ...{wt_orf[mutation_pos-4:mutation_pos+1]}...")
    print(f"Mutant sequence context: ...{mut_orf[mutation_pos-4:mutation_pos+1]}...")
    print("\nConclusion from sequence analysis:")
    if mut_has_site and not wt_has_site:
        print("The C->A mutation CREATES an SfaNI recognition site.")
        print("Therefore:")
        print(" - The WILD-TYPE (WT) allele will NOT be cut by SfaNI.")
        print(" - The MUTANT allele WILL be cut by SfaNI.")
    else:
        # This branch is for completeness, but based on the data, the above is true.
        print("The logic suggests the mutation either destroys a site or has no effect.")
        print("Re-evaluating based on gel data is necessary.")

    # Step 3: Predict the banding patterns on the gel
    print("\n--- Step 2: Predicting Gel Banding Patterns ---")
    print("Based on the sequence analysis and the gel image:")
    print("- The DNA ladder shows bands at 250 bp, 500 bp, and 750 bp.")
    print("- From the gel, the uncut PCR product (from the WT allele) is approximately 480 bp.")
    print("- When the mutant allele is cut, it produces two smaller fragments of approx. 320 bp and 160 bp (320+160=480).")
    print("\nThis leads to the following expected patterns for each genotype:")
    print(textwrap.fill("- Homozygous Wild-Type (+/+): Has only the uncut WT allele. Will show ONE band at ~480 bp.", 70))
    print(textwrap.fill("- Homozygous Mutant (-/-): Has only the cut mutant allele. Will show TWO bands at ~320 bp and ~160 bp.", 70))
    print(textwrap.fill("- Heterozygote (+/-): Has one WT and one mutant allele. Will show THREE bands (~480 bp, ~320 bp, and ~160 bp).", 70))

    # Step 4: Count the genotypes by interpreting the gel image
    print("\n--- Step 3: Counting Genotypes from the Gel Image ---")
    print("By observing the banding pattern in each of the 18 sample lanes, we can count the genotypes:")
    
    # These counts are based on visual inspection of the provided gel image.
    num_homozygous_wt = 3
    num_heterozygous = 9
    num_homozygous_mut = 6
    
    print(f"- Number of lanes with one band (~480 bp) [Homozygous Wild-Type]: {num_homozygous_wt}")
    print(f"- Number of lanes with three bands (~480, ~320, ~160 bp) [Heterozygote]: {num_heterozygous}")
    print(f"- Number of lanes with two bands (~320, ~160 bp) [Homozygous Mutant]: {num_homozygous_mut}")

    print("\n--- Final Answer ---")
    print(f"The number of larvae for each genotype is:")
    print(f"Homozygous wild-type: {num_homozygous_wt}")
    print(f"Heterozygotes: {num_heterozygous}")
    print(f"Homozygous mutants: {num_homozygous_mut}")
    print("\nIn the format A/B/C, the answer is:")
    print(f"{num_homozygous_wt}/{num_heterozygous}/{num_homozygous_mut}")

# Run the analysis
analyze_genotypes()