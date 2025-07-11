def analyze_gene_flow_consequences():
    """
    This script explains the relationship between gene flow and various population
    genetic metrics to determine which scenario is inconsistent with gene flow
    across a hybrid zone.
    """
    print("Analyzing the effect of gene flow across a hybrid zone:")
    print("========================================================\n")

    print("Core Concept: Gene flow is the movement of genes between populations.")
    print("This process makes populations more genetically similar and counteracts differentiation.\n")

    print("--- Analyzing the Options ---\n")

    # Option A: Fst
    print("A. High Fst between populations:")
    print("   - Fst measures population differentiation. A high Fst (close to 1) means populations are very different.")
    print("   - Gene flow's primary role is to reduce differentiation.")
    print("   - CONCLUSION: High gene flow and high Fst are mutually exclusive. If gene flow is occurring, Fst must be lowered. Therefore, this is the scenario that CANNOT occur.\n")

    # Option B: Dxy
    print("B. High Dxy between populations:")
    print("   - Dxy measures the average number of DNA differences between two populations, reflecting their divergence time.")
    print("   - Gene flow can occur between two populations that diverged long ago and thus have a high Dxy.")
    print("   - CONCLUSION: A high Dxy is compatible with gene flow.\n")

    # Option C: Fis
    print("C. High Fis within a population:")
    print("   - Fis measures the deficit of heterozygotes within a population (inbreeding).")
    print("   - In a hybrid zone, the mixing of two distinct groups (the 'Wahlund effect') can create a statistical deficit of heterozygotes.")
    print("   - CONCLUSION: A high Fis is compatible with gene flow in a hybrid zone.\n")
    
    # Option D: High u (mutation rate)
    print("D. High u (mutation rate) within a population:")
    print("   - The mutation rate (u) is a fundamental biological property of an organism's genome.")
    print("   - Gene flow moves alleles; it does not influence the rate at which new mutations arise.")
    print("   - CONCLUSION: A high mutation rate is compatible with gene flow.\n")

    # Option E: High Pi
    print("E. High Pi within a population:")
    print("   - Pi (nucleotide diversity) is the average genetic variation *within* a population.")
    print("   - Gene flow introduces new alleles from another population, which increases genetic variation.")
    print("   - CONCLUSION: A high Pi is not only compatible but is expected in a hybrid zone with gene flow.\n")

    print("========================================================\n")
    print("Final Determination: The presence of significant gene flow actively works to reduce Fst. Therefore, a high Fst between populations cannot occur simultaneously with high gene flow.")

analyze_gene_flow_consequences()
<<<A>>>