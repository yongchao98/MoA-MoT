def analyze_hybrid_zone_genetics():
    """
    This script analyzes population genetic metrics in the context of a hybrid zone
    with active gene flow to determine which scenario is impossible.
    """

    print("Analyzing the core question: When gene flow occurs across a hybrid zone, which of the following cannot occur?")
    print("The key concept is that gene flow makes populations more genetically similar.")
    print("-" * 50)

    # --- Analysis of Option A ---
    print("Option A: High Fst between populations")
    print("  - What is Fst? Fst is the Fixation Index, which measures genetic differentiation between populations.")
    print("  - An Fst value close to 1 indicates high differentiation (populations are genetically distinct).")
    print("  - An Fst value close to 0 indicates low differentiation (populations are genetically similar).")
    print("  - How does gene flow affect Fst? Gene flow's primary effect is to homogenize populations, which directly REDUCES Fst.")
    print("  - Conclusion: The presence of significant gene flow is fundamentally incompatible with high Fst. Therefore, high Fst CANNOT occur.")
    print("-" * 50)

    # --- Analysis of Option B ---
    print("Option B: High Dxy between populations")
    print("  - What is Dxy? Dxy measures the average absolute number of nucleotide differences between two populations.")
    print("  - It reflects the historical divergence time. If populations were separated for a long time, Dxy will be high.")
    print("  - Conclusion: Even with current gene flow, a high Dxy can persist as a 'fossil' of deep historical separation. High Dxy CAN occur.")
    print("-" * 50)

    # --- Analysis of Option C ---
    print("Option C: High Fis within a population")
    print("  - What is Fis? Fis is the Inbreeding Coefficient, measuring a deficit of heterozygotes within a population.")
    print("  - In a hybrid zone, assortative mating (like-with-like) or selection against hybrids (who are heterozygous) can cause a heterozygote deficit.")
    print("  - Conclusion: A high Fis is a common feature of hybrid zones. High Fis CAN occur.")
    print("-" * 50)

    # --- Analysis of Option D ---
    print("Option D: High u (mutation rate) within a population")
    print("  - What is u? u is the intrinsic rate at which new mutations arise.")
    print("  - Gene flow (movement of existing alleles) and mutation (creation of new alleles) are independent processes.")
    print("  - Conclusion: A population can have a high mutation rate regardless of gene flow. High u CAN occur.")
    print("-" * 50)

    # --- Analysis of Option E ---
    print("Option E: High Pi within a population")
    print("  - What is Pi (π)? Pi is nucleotide diversity, a measure of genetic variation WITHIN a population.")
    print("  - Gene flow from a different population introduces new alleles, increasing the total variation within the hybrid zone.")
    print("  - Conclusion: Mixing two divergent gene pools often leads to high Pi. High Pi CAN occur.")
    print("-" * 50)

# Run the analysis function
analyze_hybrid_zone_genetics()
<<<A>>>