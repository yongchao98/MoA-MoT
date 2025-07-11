def analyze_heritability():
    """
    Analyzes the relationship between broad-sense heritability, narrow-sense
    heritability, and the predictive power of a polygenic score (PGS).
    """

    # --- Step 1: Define knowns from the problem ---
    H2 = 0.5  # Broad-sense heritability
    # For demonstration, let's assume a total phenotypic variance (Vp)
    Vp = 100.0

    # --- Step 2: Calculate total genetic variance (Vg) ---
    # H2 = Vg / Vp  =>  Vg = H2 * Vp
    Vg = H2 * Vp

    print("--- Core Principles ---")
    print(f"Given Broad-Sense Heritability (H2) = {H2}")
    print(f"Assuming a Total Phenotypic Variance (Vp) = {Vp}")
    print(f"The Total Genetic Variance (Vg) is H2 * Vp = {H2} * {Vp} = {Vg}")
    print("\nGenetic Variance (Vg) is the sum of Additive (Va), Dominance (Vd), and Epistatic (Vi) variances.")
    print("Vg = Va + Vd + Vi")
    print("\nA standard Polygenic Score (PGS) from GWAS explains Additive Variance (Va).")
    print("Its maximum explanatory power is the Narrow-Sense Heritability (h2) = Va / Vp.")

    # --- Step 3: Demonstrate the core inequality ---
    print("\n--- The Key Relationship ---")
    print("Since Vd and Vi cannot be negative (they are variances), Va must be less than or equal to Vg.")
    print("Va <= Vg")
    print("Dividing by Vp, we get: (Va / Vp) <= (Vg / Vp)")
    print("Therefore, it is a necessary truth that: h2 <= H2")
    print(f"Since H2 = {H2}, the maximum possible value for h2 is {H2}.")
    print("This means a PGS can, at most, explain 50% of the variance.")

    # --- Step 4: Illustrate with scenarios ---
    print("\n--- Illustrative Scenarios ---")

    # Scenario 1: Non-additive effects exist (Vd > 0 or Vi > 0)
    Vd_1 = 15.0
    Vi_1 = 5.0
    Va_1 = Vg - Vd_1 - Vi_1
    h2_1 = Va_1 / Vp
    print("\nScenario A: Some genetic variance is non-additive.")
    print(f"  If Vg = {Vg}, and we assume Vd = {Vd_1} and Vi = {Vi_1}:")
    print(f"  Then Va = Vg - Vd - Vi = {Vg} - {Vd_1} - {Vi_1} = {Va_1}")
    print(f"  The resulting h2 = Va / Vp = {Va_1} / {Vp} = {h2_1}")
    print(f"  In this case, the PGS explains {h2_1:.2%}, which is less than 50%.")

    # Scenario 2: All genetic variance is additive (Vd = 0, Vi = 0)
    Vd_2 = 0.0
    Vi_2 = 0.0
    Va_2 = Vg - Vd_2 - Vi_2
    h2_2 = Va_2 / Vp
    print("\nScenario B: All genetic variance is additive (a theoretical case).")
    print(f"  If Vg = {Vg}, and we assume Vd = {Vd_2} and Vi = {Vi_2}:")
    print(f"  Then Va = Vg - Vd - Vi = {Vg} - {Vd_2} - {Vi_2} = {Va_2}")
    print(f"  The resulting h2 = Va / Vp = {Va_2} / {Vp} = {h2_2}")
    print(f"  In this case, the PGS explains {h2_2:.2%}, which is exactly 50%.")

    # --- Step 5: Final conclusion ---
    print("\n--- Final Conclusion ---")
    print("As shown, the variance explained by the PGS (h2) can be equal to 50% or less than 50%.")
    print("However, it can NEVER be MORE than 50%.")
    print("Therefore, the only statement that is necessarily true is:")
    print("A. The polygenic score can not explain more than 50% of the variance in the phenotype")

if __name__ == '__main__':
    analyze_heritability()