import sys

def heritability_analysis():
    """
    Analyzes and explains the relationship between heritability and polygenic scores.
    """
    # --- Step 1: Define parameters from the problem ---
    # We can set the total phenotypic variance (Vp) to 1.0 without loss of generality.
    Vp = 1.0
    # The problem states that broad-sense heritability (H^2) is 0.5.
    H2 = 0.5
    # From H^2 = Vg / Vp, we can calculate the total genetic variance (Vg).
    Vg = H2 * Vp

    print("--- Problem Setup ---")
    print(f"Total Phenotypic Variance (Vp) is set to: {Vp}")
    print(f"Broad-Sense Heritability (H²) is given as: {H2}")
    print(f"This implies Total Genetic Variance (Vg) = H² * Vp = {H2} * {Vp} = {Vg}")
    print("-" * 40)
    print("\nExplanation:")
    print("Total Genetic Variance (Vg) is composed of Additive (Va), Dominance (Vd), and Epistatic (Vi) variance.")
    print("Vg = Va + Vd + Vi")
    print("Narrow-sense heritability (h²) depends only on additive variance: h² = Va / Vp.")
    print("A standard polygenic score (PGS) is built on additive effects, so the maximum variance it can explain is h².")
    print("-" * 40)


    # --- Step 2: Model different scenarios for the composition of Vg ---

    # Scenario A: All genetic variance is purely additive.
    Va_scenario_A = Vg
    V_non_additive_A = Vg - Va_scenario_A
    h2_A = Va_scenario_A / Vp
    max_pgs_r2_A = h2_A

    print("\n--- Scenario A: All Genetic Variance is Additive ---")
    print("In this case, Va = Vg and non-additive variance (Vd + Vi) = 0.")
    print(f"Additive Variance (Va) = {Va_scenario_A}")
    print(f"Non-Additive Variance (Vd + Vi) = {V_non_additive_A}")
    print(f"Narrow-Sense Heritability (h²) = Va / Vp = {Va_scenario_A} / {Vp} = {h2_A}")
    print(f"Max variance a PGS can explain = h² = {max_pgs_r2_A}")
    print(f"Result: The PGS can explain up to 50% of the variance.")

    # Scenario B: A more typical case with non-additive genetic effects.
    # Let's assume Va is 70% of Vg.
    Va_scenario_B = Vg * 0.7
    V_non_additive_B = Vg - Va_scenario_B
    h2_B = Va_scenario_B / Vp
    max_pgs_r2_B = h2_B

    print("\n--- Scenario B: Non-Additive Variance is Present ---")
    print("In this case, Va < Vg. Let's assume Va is 70% of Vg.")
    print(f"Additive Variance (Va) = {Va_scenario_B:.2f}")
    print(f"Non-Additive Variance (Vd + Vi) = {V_non_additive_B:.2f}")
    print(f"Narrow-Sense Heritability (h²) = Va / Vp = {Va_scenario_B:.2f} / {Vp} = {h2_B:.2f}")
    print(f"Max variance a PGS can explain = h² = {max_pgs_r2_B:.2f}")
    print(f"Result: The PGS can explain up to 35% of the variance.")


    # --- Step 3: Draw the final conclusion ---
    print("\n" + "-" * 40)
    print("Conclusion:")
    print("The key relationship is: Variance_explained_by_PGS <= h² <= H².")
    print(f"Since H² = 0.5, it is always true that the variance explained by the PGS must be <= 0.5.")
    print("Therefore, the polygenic score can NOT explain MORE than 50% of the variance.")
    print("\nStatement A is the only one that is necessarily true.")


if __name__ == '__main__':
    heritability_analysis()
