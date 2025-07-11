def analyze_heritability_statements():
    """
    Analyzes statements about heritability and polygenic scores
    by printing a step-by-step logical explanation.
    """

    # --- Step 1: Define the core concepts from the problem ---
    print("--- Step 1: Defining Key Concepts ---")
    H_squared = 0.5
    print("Phenotypic Variance (Vp) = Genetic Variance (Vg) + Environmental Variance (Ve)")
    print(f"We are given a broad-sense heritability (H^2) of {H_squared}.")
    print("Broad-sense heritability is the proportion of phenotypic variance explained by ALL genetic factors.")
    print(f"Equation: H^2 = Vg / Vp = {H_squared}")
    print("This means that, at most, 50% of the variance in the phenotype is due to genetic factors of any kind (additive, dominance, interactions).\n")

    # --- Step 2: Establish the fundamental limit of a Polygenic Score ---
    print("--- Step 2: The Fundamental Limit for a Genetic Predictor ---")
    print("A Polygenic Score (PGS) is a predictor built using an individual's genetic data.")
    print("Therefore, the variance a PGS can explain (R^2_PGS) is fundamentally limited by the total genetic variance (Vg).")
    print("No genetic predictor, no matter how perfect, can explain more variance than is caused by genetics in the first place.")
    print("This gives us a core inequality:")
    print("R^2_PGS <= H^2")
    print(f"R^2_PGS <= {H_squared}\n")

    # --- Step 3: Evaluate each answer choice based on this logic ---
    print("--- Step 3: Evaluating The Answer Choices ---\n")

    # Analysis of Statement A
    print("A. The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print("   - This statement is a direct consequence of the fundamental limit R^2_PGS <= H^2.")
    print(f"   - Since H^2 = {H_squared}, the polygenic score cannot explain more than {H_squared * 100}% of the variance.")
    print("   - This statement is NECESSARILY TRUE.\n")

    # Analysis of Statement B
    print("B. Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print("   - A standard PGS from a GWAS captures *additive* genetic variance (Va), and its maximum explanatory power is the narrow-sense heritability (h^2 = Va / Vp).")
    print(f"   - This statement implies that h^2 would approach H^2 ({H_squared}). This only happens if all genetic variance is additive.")
    print("   - If non-additive effects (dominance, epistasis) exist, then h^2 < H^2. For example, h^2 could be 0.3, and the PGS would never reach 0.5.")
    print("   - We cannot assume non-additive effects are zero. Thus, this is NOT necessarily true.\n")

    # Analysis of Statement C
    print("C. Given an arbitrarily large GWAS, the PGS ... will not approach a variance explained of 50% due to ... non-linear effects.")
    print("   - This assumes that non-additive (non-linear) effects MUST exist.")
    print("   - It is theoretically possible, although perhaps unlikely, that all genetic variance for this trait is additive. In that case, h^2 = H^2 = 0.5.")
    print("   - Because we cannot rule out this possibility, we cannot say this statement is NECESSARILY true.\n")

    # Analysis of Statement D
    print("D. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print(f"   - The initial H^2 of {H_squared} already partitions variance into genetic (50%) and non-genetic (50%).")
    print("   - Epigenetic effects are typically modeled as part of the non-genetic (environmental) variance.")
    print(f"   - Their existence doesn't change the fact that h^2 <= H^2 ({H_squared}). It's still possible for h^2 to equal 0.5 if all genetic variance is additive.")
    print("   - This is NOT necessarily true.\n")

# Run the analysis
analyze_heritability_statements()