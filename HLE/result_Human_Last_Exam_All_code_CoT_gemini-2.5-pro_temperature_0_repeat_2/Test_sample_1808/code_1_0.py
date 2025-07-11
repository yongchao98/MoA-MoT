def analyze_gene_flow_effects():
    """
    Analyzes the effect of gene flow on various population genetic metrics
    to determine which is incompatible with it.
    """
    print("Analyzing the relationship between gene flow and population genetic metrics:\n")

    # --- Definitions and Analysis ---
    print("A. High Fst (Fixation Index):")
    print("   - Fst measures genetic differentiation between populations.")
    print("   - Gene flow homogenizes populations, moving alleles between them.")
    print("   - This action directly REDUCES Fst.")
    print("   - Therefore, high gene flow is incompatible with a high Fst. A high Fst implies low gene flow.\n")

    print("B. High Dxy (Absolute Divergence):")
    print("   - Dxy measures the average number of nucleotide differences, reflecting historical divergence time.")
    print("   - Even with current gene flow, Dxy can remain high if the populations were separated for a long time before contact.")
    print("   - Therefore, high Dxy CAN occur with gene flow.\n")

    print("C. High Fis (Inbreeding Coefficient):")
    print("   - Fis measures heterozygote deficit within a population.")
    print("   - A hybrid zone is a mix of two distinct gene pools. This mixture (the Wahlund effect) causes a statistical deficit of heterozygotes.")
    print("   - Therefore, a high Fis CAN occur with gene flow in a hybrid zone.\n")

    print("D. High u (Mutation Rate):")
    print("   - 'u' is the rate of new mutations. This is a fundamental biological process.")
    print("   - Gene flow is the movement of existing alleles and does not affect the mutation rate.")
    print("   - Therefore, a high mutation rate CAN occur with gene flow (it is an independent factor).\n")

    print("E. High Pi (Nucleotide Diversity):")
    print("   - Pi measures genetic diversity within a population.")
    print("   - Gene flow brings new alleles from one population to another, increasing the total genetic variation in the hybrid zone.")
    print("   - Therefore, a high Pi is an EXPECTED outcome of gene flow in a hybrid zone.\n")

    print("--- Conclusion ---")
    print("The only metric that is directly counteracted by gene flow is Fst.")
    print("High gene flow leads to low Fst. Thus, the occurrence of gene flow is incompatible with maintaining a high Fst.")

analyze_gene_flow_effects()