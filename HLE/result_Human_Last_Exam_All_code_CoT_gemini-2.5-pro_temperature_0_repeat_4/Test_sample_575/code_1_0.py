def analyze_heritability_scenarios():
    """
    This function models the components of genetic variance to evaluate the given statements.
    It demonstrates that the variance explained by a Polygenic Score (PGS) is always
    less than or equal to the broad-sense heritability (H^2).
    """
    # --- Constants from the problem ---
    # Broad-sense heritability (H^2) is the total proportion of variance due to genetics.
    H2_given = 0.5
    # We can assume a total phenotypic variance (Vp) of 100 for easy interpretation of percentages.
    Vp_assumed = 100.0

    print("Problem Setup:")
    print(f"Broad-Sense Heritability (H^2) = {H2_given}")
    print(f"Assumed Total Phenotypic Variance (Vp) = {Vp_assumed}")
    print("-" * 40)

    # From the setup, we can calculate the total genetic variance (Vg).
    # Vg = H^2 * Vp
    Vg = H2_given * Vp_assumed
    print(f"Calculated Total Genetic Variance (Vg) = {H2_given} * {Vp_assumed} = {Vg}")
    print("-" * 40)

    # --- Scenario 1: All genetic variance is additive ---
    # This scenario tests if statement C ("will not approach 50%") is necessarily true.
    print("\nScenario 1: Assume all genetic variance is additive.")
    fraction_additive_1 = 1.0
    Va_1 = Vg * fraction_additive_1
    V_non_additive_1 = Vg - Va_1
    h2_1 = Va_1 / Vp_assumed
    R2_PGS_1 = h2_1

    print(f"Additive Variance (Va) = {Vg} * {fraction_additive_1} = {Va_1}")
    print(f"Non-Additive Variance (Vd + Vi) = {Vg} - {Va_1} = {V_non_additive_1}")
    print(f"Narrow-Sense Heritability (h^2 = Va / Vp) = {Va_1} / {Vp_assumed} = {h2_1}")
    print(f"In this case, an ideal PGS would explain {R2_PGS_1*100:.0f}% of the variance.")
    print("Since it's possible for the PGS to approach 50%, statement C is not necessarily true.")
    print("-" * 40)

    # --- Scenario 2: Some genetic variance is non-additive ---
    # This scenario tests if statement B ("will approach 50%") is necessarily true.
    print("\nScenario 2: Assume 60% of genetic variance is additive.")
    fraction_additive_2 = 0.6
    Va_2 = Vg * fraction_additive_2
    V_non_additive_2 = Vg - Va_2
    h2_2 = Va_2 / Vp_assumed
    R2_PGS_2 = h2_2

    print(f"Additive Variance (Va) = {Vg} * {fraction_additive_2} = {Va_2}")
    print(f"Non-Additive Variance (Vd + Vi) = {Vg} - {Va_2} = {V_non_additive_2}")
    print(f"Narrow-Sense Heritability (h^2 = Va / Vp) = {Va_2} / {Vp_assumed} = {h2_2}")
    print(f"In this case, an ideal PGS would explain {R2_PGS_2*100:.0f}% of the variance.")
    print("Since the PGS explains less than 50%, statement B is not necessarily true.")
    print("-" * 40)

    # --- Final Conclusion ---
    print("\nFinal Conclusion:")
    print(f"In Scenario 1, PGS variance explained = {R2_PGS_1}. Is this <= H^2 ({H2_given})? {R2_PGS_1 <= H2_given}")
    print(f"In Scenario 2, PGS variance explained = {R2_PGS_2}. Is this <= H^2 ({H2_given})? {R2_PGS_2 <= H2_given}")
    print("\nIn all possible scenarios, the variance explained by a PGS is capped by the total genetic variance.")
    print(f"The total genetic variance is given as {H2_given*100:.0f}% of the total phenotypic variance (H^2 = {H2_given}).")
    print("Therefore, a polygenic score cannot explain more than 50% of the variance.")
    print("Statement A is the only one that is necessarily true.")

if __name__ == '__main__':
    analyze_heritability_scenarios()