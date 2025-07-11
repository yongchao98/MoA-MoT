def explain_hybrid_zone_genetics():
    """
    Explains why high Fst cannot occur with significant gene flow in a hybrid zone.
    This function will print the reasoning for each option and then state the final conclusion.
    """

    print("Analyzing the effects of gene flow in a hybrid zone on various population genetic statistics:\n")

    # A. High Fst
    print("A. High Fst (Population Differentiation):")
    print("   - Fst measures how genetically different two populations are.")
    print("   - Gene flow is the movement of genes between populations, which makes them more similar.")
    print("   - Therefore, gene flow actively REDUCES Fst.")
    print("   - If significant gene flow is occurring, it is impossible for populations to remain highly differentiated.")
    print("   - CONCLUSION: High Fst and high gene flow are mutually exclusive. This CANNOT occur.\n")

    # B. High Dxy
    print("B. High Dxy (Absolute Divergence):")
    print("   - Dxy measures the average number of DNA differences, reflecting the time since two populations split.")
    print("   - If populations were isolated for a long time before forming the hybrid zone, Dxy will be high.")
    print("   - Gene flow will slowly decrease Dxy, but it doesn't erase the historical signal immediately.")
    print("   - CONCLUSION: High Dxy CAN occur, as it reflects deep history.\n")

    # C. High Fis
    print("C. High Fis (Inbreeding Coefficient):")
    print("   - High Fis indicates a deficit of heterozygotes.")
    print("   - In a hybrid zone, you have a mix of parent population 1 (e.g., AA), parent population 2 (e.g., aa), and hybrids (Aa).")
    print("   - This mixture, known as the Wahlund Effect, creates a statistical deficit of heterozygotes in the combined sample.")
    print("   - CONCLUSION: High Fis is EXPECTED to occur in a hybrid zone.\n")

    # D. High u (mutation rate)
    print("D. High u (mutation rate):")
    print("   - The mutation rate (u or µ) is a basic biological rate at which DNA changes.")
    print("   - It is not directly affected by the movement of individuals between populations (gene flow).")
    print("   - A population could have a high mutation rate for other reasons (e.g., exposure to mutagens).")
    print("   - CONCLUSION: High µ CAN occur, as it is independent of gene flow.\n")
    
    # E. High Pi
    print("E. High Pi (Nucleotide Diversity):")
    print("   - Pi measures genetic variation WITHIN a population.")
    print("   - Gene flow from a different population introduces new alleles.")
    print("   - This mixing of two diverse gene pools increases the overall variation in the hybrid zone.")
    print("   - CONCLUSION: High Pi is EXPECTED to occur.\n")

    print("--------------------------------------------------")
    print("Final Conclusion: The one condition that is directly counteracted by gene flow is high Fst.")
    print("--------------------------------------------------")

# Execute the explanation function
explain_hybrid_zone_genetics()