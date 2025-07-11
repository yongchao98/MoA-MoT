def solve_heritability_problem():
    """
    This script logically evaluates statements about heritability and polygenic scores.
    """
    
    # --- Step 1 & 2: Define terms and relationships ---
    H_squared = 0.5
    
    print("### Analysis of the Genetics Problem ###\n")
    print("Let's break down the concepts based on the provided information.")
    
    print("1. Given Information:")
    print(f"   - Broad-sense heritability (H²), the total proportion of phenotypic variance due to all genetic factors, is {H_squared}.")
    
    print("\n2. Key Definitions:")
    print("   - Total Phenotypic Variance (Vp) = Genetic Variance (Vg) + Environmental Variance (Ve)")
    print("   - Genetic Variance (Vg) = Additive (Va) + Dominance (Vd) + Epistatic (Vi) Variance")
    print(f"   - Broad-Sense Heritability (H²) = Vg / Vp = {H_squared}")
    print("   - Narrow-Sense Heritability (h²) = Va / Vp")

    print("\n3. Polygenic Score (PGS) Limitation:")
    print("   - A standard PGS from a GWAS primarily captures linear, additive genetic effects.")
    print("   - Therefore, the maximum possible variance a PGS can explain is equal to the narrow-sense heritability (h²).\n")

    # --- Step 3: The Core Inequality ---
    print("--- The Core Logical Deduction ---")
    print("By definition, the Additive Variance (Va) is a component of the total Genetic Variance (Vg).")
    print("This means Va must be less than or equal to Vg (Va <= Vg), since dominance (Vd) and epistatic (Vi) variance cannot be negative.")
    print("If we divide both sides by the total Phenotypic Variance (Vp), the relationship holds:")
    print("=> (Va / Vp) <= (Vg / Vp)")
    print("This translates directly to our heritability measures:")
    print("=> h² <= H²")
    print("Therefore, the maximum variance a PGS can explain (which is h²) must be less than or equal to the broad-sense heritability, H².\n")

    # --- Step 4: Evaluate the statements ---
    print("--- Evaluating the Answer Choices ---\n")

    print("Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print(f"   - Reasoning: The maximum variance a PGS can explain is h². We know h² <= H², and H² = {H_squared}.")
    print(f"   - The final 'equation' or inequality is: PGS_explained_variance <= h² <= {H_squared}")
    print("   - This statement holds true no matter the value of h². It correctly states the upper bound.")
    print("   - Verdict: NECESSARILY TRUE.\n")

    print("Statement B: Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print(f"   - Reasoning: This would only be true if h² approaches H² (i.e., h² = {H_squared}). This requires non-additive genetic variance (Vd + Vi) to be zero.")
    print("   - The problem does not guarantee this condition.")
    print("   - Verdict: NOT necessarily true.\n")

    print("Statement C: Given an arbitrarily large GWAS, the PGS ... will not approach a variance explained of 50% due to gene-gene interactions and other non-linear effects...")
    print(f"   - Reasoning: This claims that h² will be strictly less than {H_squared}. This requires non-additive variance (Vd + Vi) to be greater than zero.")
    print("   - While likely, the problem does not guarantee the existence of non-additive effects.")
    print("   - Verdict: NOT necessarily true.\n")

    print("Statement D: The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print(f"   - Reasoning: Epigenetic variance is part of the total phenotypic variance (Vp). The value H² = {H_squared} already accounts for this.")
    print(f"   - It is possible that h² = H² = {H_squared}, with all non-genetic variance being environmental or epigenetic.")
    print("   - Verdict: NOT necessarily true.\n")


# Execute the analysis
solve_heritability_problem()