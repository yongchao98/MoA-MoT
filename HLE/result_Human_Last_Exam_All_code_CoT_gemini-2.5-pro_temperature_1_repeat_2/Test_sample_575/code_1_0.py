import math

def analyze_heritability():
    """
    Analyzes the maximum variance explained by a Polygenic Score (PGS)
    based on a given broad-sense heritability (H^2).
    """
    # --- Given Parameters ---
    # Broad-sense heritability (H^2) is the total proportion of variance due to genetics.
    H_squared = 0.5

    # To make calculations concrete, let's assume a total phenotypic variance (Vp).
    # The specific value doesn't matter, as heritability is a ratio. We'll use 100.
    V_p = 100.0
    print("--- Basic Definitions ---")
    print(f"Assumed Total Phenotypic Variance (Vp): {V_p}")
    print(f"Given Broad-Sense Heritability (H^2): {H_squared}")

    # --- Fundamental Limit ---
    # Total Genetic Variance (Vg) is Vp * H^2.
    V_g = V_p * H_squared
    print("\n--- Fundamental Limit ---")
    print(f"The total variance due to all genetic factors (Vg) is calculated as:")
    print(f"Vg = Vp * H^2 = {V_p} * {H_squared} = {V_g}")
    print("\nThe broad-sense heritability (H^2) represents the ABSOLUTE MAXIMUM variance")
    print("that can be explained by any predictor based on genetics.")
    print(f"Maximum possible R^2 from genetics = Vg / Vp = {V_g} / {V_p} = {V_g / V_p}")
    print("\nThis directly means that a polygenic score cannot explain more than 50% of the variance.")
    print("-" * 40)

    # --- Illustrative Scenarios ---
    # A standard PGS built from GWAS primarily captures additive genetic effects.
    # Its performance is limited by narrow-sense heritability (h^2 = Va / Vp).
    # Let's explore two possibilities for the genetic architecture.

    print("\n--- Scenario 1: All genetic variance is additive ---")
    # In this case, dominance (Vd) and epistasis (Vi) are zero.
    # Va = Vg, Vd = 0, Vi = 0
    V_a_1 = V_g
    h_squared_1 = V_a_1 / V_p
    print("Assume all genetic variance is additive (Va = Vg).")
    print(f"Additive Genetic Variance (Va) = {V_a_1}")
    print(f"Narrow-Sense Heritability (h^2) is calculated as:")
    print(f"h^2 = Va / Vp = {V_a_1} / {V_p} = {h_squared_1}")
    print(f"In this case, the theoretical max R^2 for a standard PGS is {h_squared_1*100}%.")

    print("\n--- Scenario 2: Genetic variance includes non-additive effects ---")
    # Assume half of the genetic variance is additive and half is non-additive.
    # Va < Vg
    V_a_2 = V_g / 2.0
    V_non_additive = V_g / 2.0
    h_squared_2 = V_a_2 / V_p
    print("Assume 50% of genetic variance is additive and 50% is non-additive.")
    print(f"Total Genetic Variance (Vg) = {V_g}")
    print(f"Additive Genetic Variance (Va) = {V_a_2}")
    print(f"Non-Additive Genetic Variance (Vd + Vi) = {V_non_additive}")
    print(f"Narrow-Sense Heritability (h^2) is calculated as:")
    print(f"h^2 = Va / Vp = {V_a_2} / {V_p} = {h_squared_2}")
    print(f"In this case, the theoretical max R^2 for a standard PGS is {h_squared_2*100}%.")
    print("-" * 40)

    # --- Final Conclusion ---
    print("\nConclusion:")
    print("In all scenarios, the variance explained by a PGS is less than or equal to H^2 (0.5).")
    print(f"Scenario 1 Limit: {h_squared_1} <= {H_squared}")
    print(f"Scenario 2 Limit: {h_squared_2} <= {H_squared}")
    print("Therefore, the statement 'The polygenic score can not explain more than 50% of the variance in the phenotype' is necessarily true.")


if __name__ == '__main__':
    analyze_heritability()