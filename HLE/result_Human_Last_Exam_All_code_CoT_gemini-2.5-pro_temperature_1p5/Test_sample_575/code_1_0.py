def explain_heritability():
    """
    Illustrates the principles of heritability and polygenic score limits.
    """

    # --- Given Information ---
    # Broad-sense heritability (H_squared) is the total proportion of phenotypic variance
    # explained by all genetic factors (additive, dominance, epistasis).
    H_squared = 0.5

    # For illustration, let's assume a total phenotypic variance (V_p) of 100 units.
    V_p = 100.0

    # --- Analysis of Statement A ---
    print("--- Analyzing Statement A ---")
    print("A. The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print("\nStep 1: Define the relationship between variances.")
    print(f"Broad-sense Heritability (H²): {H_squared}")
    print(f"Total Phenotypic Variance (Vp): {V_p} (assumed for illustration)")

    # Calculate the total genetic variance (V_g).
    # H² = Vg / Vp  =>  Vg = H² * Vp
    V_g = H_squared * V_p
    print("\nStep 2: Calculate the total genetic variance (Vg).")
    print(f"Vg = H² * Vp")
    print(f"Vg = {H_squared} * {V_p} = {V_g}")
    print("\nThis means that all genetic factors combined account for 50.0 units of the total 100.0 units of variance.")

    print("\nStep 3: State the fundamental limit.")
    print("A polygenic score (PGS) predicts a phenotype using only genetic information.")
    print("Therefore, the variance it can explain (let's call it R²_PGS) is limited by the total genetic variance (Vg).")
    print("The maximum possible variance a PGS can explain is Vg.")
    print(f"Expressed as a proportion of total variance: R²_PGS <= Vg / Vp")
    print(f"Since H² = Vg / Vp, it means: R²_PGS <= H²")
    print(f"R²_PGS <= {H_squared}")
    print("\nConclusion: The PGS can not explain more than 50% of the variance. Statement A is necessarily TRUE.\n")

    # --- Analysis of Statements B and C ---
    print("--- Analyzing Statements B & C ---")
    print("B/C concern a standard linear PGS approaching or not approaching 50% variance explained.")
    # Genetic variance (Vg) is composed of V_a (additive), V_d (dominance), and V_i (epistasis).
    # Vg = V_a + V_d + V_i
    # We know Vg = 50.0. Let's consider a plausible scenario where non-additive effects exist.
    V_a_scenario = 40.0
    V_d_scenario = 8.0
    V_i_scenario = 2.0
    Vg_check = V_a_scenario + V_d_scenario + V_i_scenario
    print(f"\nLet's assume a plausible breakdown of Vg ({Vg_check}):")
    print(f"  - Additive Variance (Va) = {V_a_scenario}")
    print(f"  - Dominance Variance (Vd) = {V_d_scenario}")
    print(f"  - Epistatic Variance (Vi) = {V_i_scenario}")

    # Narrow-sense heritability (h²) is the proportion of variance due to ADDITIVE effects.
    # A standard linear PGS's maximum explanatory power approaches h².
    h_squared = V_a_scenario / V_p
    print("\nNarrow-sense heritability (h²) = Va / Vp")
    print(f"h² = {V_a_scenario} / {V_p} = {h_squared}")
    print(f"In this scenario, a standard PGS would approach explaining {h_squared*100}% of the variance, not {H_squared*100}%.")
    print("Because it's possible that h² < H², we cannot say statement B is necessarily true.")
    print("Conversely, because it's possible that Vd and Vi are 0 (making h² = H²), we cannot say statement C is necessarily true.")

# Run the explanation
explain_heritability()