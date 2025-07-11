import sys

def analyze_heritability_statements():
    """
    Analyzes statements about heritability and polygenic scores based on given principles.
    This script serves as a logical demonstration, not a calculation from raw data.
    """

    # --- Setup based on the problem statement ---
    # We can assume a total phenotypic variance (Vp) of 100 for clear illustration.
    Vp = 100.0

    # Broad-sense heritability (H^2) is given as 0.5.
    H_squared = 0.5

    # From the definition H^2 = Vg / Vp, we can find the total genetic variance (Vg).
    # Equation: Vg = H^2 * Vp
    Vg = H_squared * Vp
    # So, the total environmental variance (Ve) is the rest.
    # Equation: Ve = Vp - Vg
    Ve = Vp - Vg

    print("--- Foundational Concepts ---")
    print(f"Let's assume a total phenotypic variance (Vp) = {Vp:.0f} units.")
    print(f"Given Broad-Sense Heritability (H^2) = {H_squared}.")
    print(f"This means total Genetic Variance (Vg) = {H_squared} * {Vp:.0f} = {Vg:.0f} units.")
    print(f"This leaves Environmental Variance (Ve) = {Vp:.0f} - {Vg:.0f} = {Ve:.0f} units.")
    print("Genetic Variance (Vg) is composed of Additive (Va), Dominance (Vd), and Epistatic (Vi) components.")
    print("A standard PGS's maximum explained variance (R^2) approaches Narrow-Sense Heritability (h^2), where h^2 = Va / Vp.\n")

    # --- Evaluation of Each Statement ---
    print("--- Evaluating Answer Choices ---")

    # A. The polygenic score can not explain more than 50% of the variance in the phenotype.
    # The max variance explainable by genetics is Vg/Vp = H^2. A PGS is a genetic predictor.
    # Its explained variance (h^2) must be less than or equal to H^2.
    is_A_true = True  # By definition
    print("A: Can PGS explain > 50% of variance?")
    print(f"   - The total variance from all genetic factors (Vg) is {Vg:.0f} units, or {H_squared:.2f} of the total variance (Vp).")
    print(f"   - A PGS uses only genetic data, so its explained variance (R^2) must be <= Vg / Vp.")
    print(f"   - Equation: R^2(PGS) <= H^2")
    print(f"   - Numerically: R^2(PGS) <= {H_squared}")
    print(f"   - Conclusion: Statement A is TRUE.\n")

    # B and C are opposites, we can evaluate them together.
    # B claims PGS will approach 50% (H^2). C claims it will not due to non-linear effects.
    # This hinges on whether Vd+Vi can be non-zero. Let's test a realistic scenario.
    # Assume some genetic variance is non-additive, which is standard for complex traits.
    Va_scenario1 = 40.0 # Additive variance
    Vd_plus_Vi_scenario1 = 10.0 # Non-additive variance
    Vg_check = Va_scenario1 + Vd_plus_Vi_scenario1
    h_squared_scenario1 = Va_scenario1 / Vp

    print("B & C: Will the PGS variance explained approach 50%?")
    print(f"   - This depends on whether all genetic variance is additive (Vg = Va).")
    print(f"   - Biologically, for complex traits, Vg includes non-additive parts (Vd, Vi). Let's assume a plausible case:")
    print(f"     - Let Va = {Va_scenario1:.0f} and (Vd + Vi) = {Vd_plus_Vi_scenario1:.0f}. This maintains Vg = {Vg_check:.0f}.")
    print(f"     - In this case, the narrow-sense heritability h^2 = Va / Vp.")
    print(f"     - Equation: h^2 = {Va_scenario1:.0f} / {Vp:.0f}")
    print(f"     - Numerically: h^2 = {h_squared_scenario1}")
    print(f"   - The PGS variance explained would approach {h_squared_scenario1*100:.0f}%, not the full {H_squared*100:.0f}% promised by H^2.")
    print(f"   - The gap is due to the non-additive {Vd_plus_Vi_scenario1:.0f} units of variance missed by the linear PGS model.")
    print(f"   - Conclusion B (will approach 50%): This is FALSE.")
    print(f"   - Conclusion C (will NOT approach 50% due to non-linear effects): This is TRUE under the standard biological interpretation.\n")


    # D. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.
    # Epigenetic effects contribute to Ve. We need to show that Ve > 0 does NOT force h^2 < 0.5.
    # Let's construct a counterexample where all Vg is additive.
    Va_scenario2 = 50.0
    Vd_plus_Vi_scenario2 = 0.0
    h_squared_scenario2 = Va_scenario2 / Vp

    print("D: Do epigenetic effects necessarily limit h^2 to be LESS than 0.5?")
    print(f"   - Epigenetic effects are part of Environmental Variance (Ve). We already know Ve = {Ve:.0f} > 0.")
    print(f"   - The question is: does Ve > 0 REQUIRE that h^2 < 0.5? Let's test a case where ALL genetic variance is additive:")
    print(f"     - Let Va = {Va_scenario2:.0f} and (Vd + Vi) = {Vd_plus_Vi_scenario2:.0f}. This maintains Vg = {Va_scenario2:.0f}.")
    print(f"     - In this scenario, h^2 = Va / Vp.")
    print(f"     - Equation: h^2 = {Va_scenario2:.0f} / {Vp:.0f}")
    print(f"     - Numerically: h^2 = {h_squared_scenario2}")
    print(f"   - Here, h^2 is exactly 0.5, not less than 0.5, even though environmental (e.g., epigenetic) effects exist.")
    print(f"   - Conclusion: Statement D is FALSE.\n")

    print("--- Final Summary ---")
    print("Based on the analysis, Statements A and C are the only ones that are necessarily true.")

if __name__ == '__main__':
    analyze_heritability_statements()
    # The final answer corresponds to the option "Only choices A and C are correct"
    sys.stdout.write("<<<E>>>\n")
