import sys

def analyze_heritability_problem():
    """
    This script analyzes the provided genetics problem statement by statement
    to determine the necessarily true conclusion.
    """
    # --- Step 1 & 2: Define terms and establish relationships ---
    print("--- Problem Setup and Definitions ---")
    
    # H² is Broad-sense heritability, given in the problem.
    H2 = 0.5
    
    # Let's assume total phenotypic variance (Vp) is 1.0 for simplicity.
    # This makes other variances directly interpretable as proportions.
    Vp = 1.0

    # From the definition H² = Vg / Vp, we can find the total genetic variance (Vg).
    Vg = H2 * Vp
    
    # The equation for total genetic variance (Vg) is:
    # Vg = Va + Vd + Vi
    # where Va = Additive variance, Vd = Dominance variance, Vi = Epistatic (interaction) variance.
    
    # Narrow-sense heritability (h²) is the proportion of variance due to *additive* effects only.
    # The equation for h² is: h² = Va / Vp
    
    # A standard Polygenic Score (PGS) is a linear model that sums additive effects.
    # The variance explained by an ideal PGS (R²_PGS) approaches h².
    
    print(f"Given Broad-Sense Heritability (H²): {H2}")
    print(f"Let's assume Total Phenotypic Variance (Vp): {Vp}")
    print(f"This implies Total Genetic Variance (Vg) = H² * Vp = {H2} * {Vp} = {Vg:.1f}")
    print("The variance explained by a PGS (R²_PGS) is limited by narrow-sense heritability (h² = Va/Vp).")
    print("And h² is limited by H², because Va <= (Va + Vd + Vi). So, R²_PGS <= h² <= H².")
    print("-" * 40)

    # --- Step 3 & 4: Analyze each statement with scenarios ---
    print("--- Analysis of Answer Choices ---")

    # A. The polygenic score can not explain more than 50% of the variance in the phenotype
    print("\n[Analysis of A]")
    print("Is it true that R²_PGS <= 0.5?")
    print(f"The absolute ceiling for any genetic prediction is the total genetic variance, which is H² = {H2} (or 50%).")
    print("Since the PGS is a genetic predictor, its performance is fundamentally capped by H².")
    print("Conclusion: Statement A is TRUE.")

    # C. Given an arbitrarily large GWAS, the PGS ... will not approach ... 50% due to non-linear effects...
    print("\n[Analysis of C]")
    print("Is it true that R²_PGS will approach a value < 0.5?")
    print("This would be true if h² < H². This happens if non-additive variance (Vd or Vi) is greater than zero.")
    # For a "polygenic" trait in a diploid organism, the existence of non-additive effects is a standard biological assumption.
    # Let's model this realistic scenario:
    Va_realistic = 0.3
    Vd_realistic = 0.15
    Vi_realistic = 0.05
    Vg_realistic_eq = f"{Va_realistic} + {Vd_realistic} + {Vi_realistic}"
    Vg_realistic_val = Va_realistic + Vd_realistic + Vi_realistic
    h2_realistic = Va_realistic / Vp
    
    print("In a realistic polygenic trait scenario:")
    print(f"  Let Va={Va_realistic}, Vd={Vd_realistic}, Vi={Vi_realistic}")
    print(f"  The total genetic variance Vg = {Vg_realistic_eq} = {Vg_realistic_val:.1f}")
    print(f"  This is consistent with the given H² = {Vg_realistic_val / Vp:.1f}.")
    print(f"  In this case, h² = Va / Vp = {Va_realistic} / {Vp} = {h2_realistic:.1f}")
    print(f"  The PGS would approach explaining {h2_realistic*100:.0f}% of variance, which is less than 50%.")
    print("Conclusion: Under standard biological assumptions for a polygenic trait, Statement C is TRUE.")

    # B. Given an arbitrarily large GWAS, the PGS will approach a variance explained of 50%
    print("\n[Analysis of B]")
    print("Is it true that R²_PGS will approach 0.5?")
    print("This requires h² = H² = 0.5. This only happens if all genetic variance is purely additive (Vd=0 and Vi=0).")
    Va_special = 0.5
    Vd_special = 0.0
    Vi_special = 0.0
    Vg_special_eq = f"{Va_special} + {Vd_special} + {Vi_special}"
    Vg_special_val = Va_special + Vd_special + Vi_special
    print(f"  This special case requires Vg = {Vg_special_eq} = {Vg_special_val:.1f}")
    print("Since this purely additive case is not guaranteed (and is biologically unlikely), this statement is not *necessarily* true.")
    print("Conclusion: Statement B is FALSE.")

    # D. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5
    print("\n[Analysis of D]")
    print("Do epigenetic effects force h² < 0.5?")
    print("Epigenetic effects contribute to the overall variance partitioning but do not change the core definitions.")
    print(f"H² is given as 0.5. We know h² <= H², so h² <= 0.5. But h² can still equal 0.5 if all genetic variance is additive.")
    print("The source of the non-heritable variance (be it environment, epigenetics, etc.) doesn't change this relationship.")
    print("Conclusion: Statement D is FALSE.")
    print("-" * 40)
    
    # --- Step 5: Synthesize and Conclude ---
    print("--- Final Conclusion ---")
    print("Statements A and C are correct.")
    print("A is true by definition. C is true under the standard interpretation of a 'polygenic trait'.")
    print("This corresponds to option E.")

analyze_heritability_problem()