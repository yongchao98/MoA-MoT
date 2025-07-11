def analyze_heritability_statement():
    """
    Analyzes the relationship between broad-sense heritability, narrow-sense
    heritability, and polygenic scores to determine the correct statement.
    """

    # --- Step 1: Define the core concepts and the given information ---
    H2 = 0.5  # Broad-sense heritability from the problem

    print("--- Analysis of Quantitative Genetics Principles ---")
    print("1. Broad-Sense Heritability (H^2):")
    print("   - H^2 = Vg / Vp (Total Genetic Variance / Total Phenotypic Variance)")
    print(f"   - The problem states H^2 = {H2}.\n")

    print("2. Components of Genetic Variance (Vg):")
    print("   - Vg = Va + Vd + Vi")
    print("     - Va: Additive genetic variance (effects of alleles sum up linearly)")
    print("     - Vd: Dominance variance (non-linear interactions between alleles at the same locus)")
    print("     - Vi: Epistatic variance (non-linear interactions between alleles at different loci)\n")

    print("3. Narrow-Sense Heritability (h^2):")
    print("   - h^2 = Va / Vp (Additive Genetic Variance / Total Phenotypic Variance)")
    print("   - This represents the proportion of variance passed from parents to offspring.\n")

    print("4. Polygenic Score (PGS):")
    print("   - A standard PGS is built from GWAS results, which measure the linear (additive) effect of genes.")
    print("   - Therefore, the variance a PGS can explain is an estimate of h^2.\n")

    # --- Step 2: Establish the key mathematical inequality ---
    print("--- The Final Equation of Relationships ---")
    print("Based on the definitions, we have a clear hierarchy:")
    print("Variance explained by PGS <= h^2 <= H^2")
    print("\nSubstituting the value from the problem:")
    print(f"Variance explained by PGS <= h^2 <= {H2}")
    print("This inequality is the key to solving the problem. The maximum possible variance a PGS can explain is 50%.\n")

    # --- Step 3: Evaluate each answer choice ---
    print("--- Evaluating the Answer Choices ---")

    print("A. The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print(f"   - Verdict: TRUE. This is a direct consequence of the inequality 'Variance explained by PGS <= {H2}'.")
    print("   - The PGS is fundamentally limited by the total genetic variance, and even more so by the additive portion of it.\n")

    print("B. Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print(f"   - Verdict: NOT NECESSARILY TRUE. The PGS will approach h^2, not H^2.")
    print(f"   - It will only approach {H2} in the specific case where all genetic variance is additive (Vd=0 and Vi=0), which is not guaranteed.\n")

    print("C. Given an arbitrarily large GWAS, the PGS ... will not approach a variance explained of 50% due to non-linear effects...")
    print(f"   - Verdict: NOT NECESSARILY TRUE. This is the opposite of B.")
    print(f"   - If it happens that all genetic variance is additive (Vd=0 and Vi=0), then h^2 equals H^2, and the PGS *would* approach {H2}.\n")

    print("D. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print(f"   - Verdict: FALSE. The value H^2={H2} is a given ratio. We already know h^2 <= H^2.")
    print(f"   - Epigenetic effects are part of what determines the overall H^2 value, but they do not add a new constraint on h^2 beyond what is already established.\n")

    print("E. None of the other answer choices are correct.")
    print("   - Verdict: FALSE, because statement A is necessarily true.\n")

# Run the analysis
analyze_heritability_statement()

<<<A>>>