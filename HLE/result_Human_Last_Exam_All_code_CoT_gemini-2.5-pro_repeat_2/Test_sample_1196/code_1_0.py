def illustrate_misleading_association():
    """
    This function demonstrates how Tag SNPs can lead to a misleading association in a GWAS.
    """
    # Step 1: Define two different haplotypes (inherited blocks of DNA).
    # One haplotype contains a variant that causes a trait, the other does not.
    # The full sequence is: [SNP1, SNP2, SNP3(Causal Variant), SNP4]
    causal_haplotype =     ['A', 'T', 'G', 'C']
    non_causal_haplotype = ['G', 'C', 'A', 'T']

    # The 'G' at the 3rd position in the causal_haplotype is the true cause of the trait.
    causal_variant = causal_haplotype[2]
    causal_variant_position = 3

    # Step 2: Define the Tag SNPs used in the study.
    # In a typical GWAS, we don't sequence the whole genome. We genotype specific
    # "Tag SNPs" that act as markers for different haplotypes.
    # Let's say our study only genotypes SNP1 and SNP4.
    tag_snp_indices = [0, 4] # Using 1-based indexing for the printout

    # Determine the alleles of the Tag SNPs for the causal haplotype.
    tags_for_causal_haplotype = [causal_haplotype[0], causal_haplotype[3]]

    # Step 3: Explain the GWAS result.
    # In our simulated study, individuals with the `causal_haplotype` get the trait.
    # Since we only measure the Tag SNPs, we will find a strong statistical association
    # between the trait and the combination of tag alleles ['A', 'C'].
    print("--- Simulation of a Misleading GWAS Association ---")
    print(f"A full haplotype associated with a trait exists: {causal_haplotype}")
    print(f"The true CAUSAL variant is '{causal_variant}' at position {causal_variant_position}.")
    print("-" * 50)
    print("However, our GWAS only measures Tag SNPs to identify haplotypes.")
    print(f"The Tag SNPs are at positions {tag_snp_indices[0]+1} and {tag_snp_indices[1]-1}.")
    print(f"For the causal haplotype, these Tag SNPs have the alleles: {tags_for_causal_haplotype}")
    print("-" * 50)
    print("GWAS Finding:")
    print(f"A strong statistical association is found between the trait and the Tag SNP combination {tags_for_causal_haplotype}.")
    print("\nConclusion:")
    print("This finding is 'misleading' because the Tag SNPs themselves do not cause the trait.")
    print(f"They are simply markers that are inherited along with the true, unobserved causal variant ('{causal_variant}').")
    print("This perfectly illustrates Choice E: 'Three Tag SNPs predicting all alleles in an inherited haplotype.'")

# Run the illustration
illustrate_misleading_association()