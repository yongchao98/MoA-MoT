def analyze_heritability():
    """
    Analyzes statements about heritability and polygenic scores based on given premises.
    """

    # --- 1. Define Premises from the Problem ---
    # We are given the broad-sense heritability (H2).
    # Let's normalize total phenotypic variance (Vp) to 1.0 for clarity.
    Vp = 1.0
    H2 = 0.5  # Broad-sense heritability: Vg / Vp = 0.5
    Vg = H2 * Vp # Total genetic variance

    # From the definition H2 = Vg/Vp, we know that the total genetic variance Vg is 0.5 * Vp.
    print("--- Problem Setup ---")
    print(f"Total Phenotypic Variance (Vp) is set to: {Vp}")
    print(f"Broad-Sense Heritability (H2 = Vg/Vp) is given as: {H2}")
    print(f"This implies Total Genetic Variance (Vg) is: {Vg}")
    print(f"And Environmental Variance (Ve) is: {Vp - Vg}")
    print("\nGenetic Variance (Vg) is composed of Additive (Va), Dominance (Vd), and Epistatic (Vi) components.")
    print("Vg = Va + Vd + Vi")
    print("A standard Polygenic Score (PGS) variance explained approaches h2 = Va/Vp in the limit.\n")

    # --- 2. Evaluate Each Statement ---
    print("--- Evaluating the Answer Choices ---")

    # Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype.
    print("\nA. The PGS cannot explain more than 50% of the phenotypic variance.")
    # The variance explained by a PGS (let's call it V_pgs) is derived from genetic information.
    # Therefore, V_pgs must be less than or equal to the total genetic variance, Vg.
    # So, the proportion of variance explained by the PGS is V_pgs / Vp <= Vg / Vp.
    # Since Vg / Vp = H2 = 0.5, the PGS can explain at most 50%.
    print(f"   - Logic: The variance explained by any genetic model (like a PGS) cannot exceed the total genetic variance (Vg).")
    print(f"   - Equation: (Variance Explained by PGS) / Vp  <=  Vg / Vp")
    print(f"   - Calculation: (Variance Explained by PGS) / {Vp} <= {Vg} / {Vp}, which simplifies to <= {H2}")
    print("   - Conclusion: This statement is NECESSARILY TRUE.\n")

    # Statement B: Given an arbitrarily large GWAS, the PGS will approach a variance explained of 50%.
    print("B. The PGS will approach explaining 50% of the variance.")
    print("   - Logic: A standard PGS's performance limit is narrow-sense heritability (h2 = Va/Vp), not H2.")
    print("   - This is only true if all genetic variance is additive (Va = Vg).")
    # Example where it's not true:
    Va_scenario = 0.4
    Vd_Vi_scenario = Vg - Va_scenario # Non-additive variance
    h2_scenario = Va_scenario / Vp
    print(f"   - Counterexample: If Vg = {Vg} but Va = {Va_scenario} (and Vd+Vi = {Vd_Vi_scenario:.1f}), then the PGS would approach h2 = {h2_scenario:.1f} (or 40%), not 50%.")
    print("   - Conclusion: This statement is NOT necessarily true.\n")

    # Statement C: Given an arbitrarily large GWAS, the PGS will not approach ... 50% due to non-linear effects...
    print("C. The PGS will not approach 50% due to non-linear effects.")
    print("   - Logic: This statement assumes non-linear genetic effects (Vd or Vi) MUST exist.")
    print("   - The problem allows for a 'theoretically ideal' population where all genetic variance could be additive.")
    # Example where the statement is false:
    Va_scenario_2 = Vg
    print(f"   - Counterexample: If all genetic variance is additive, Va = Vg = {Vg}. Then h2 = H2 = {H2}, and the PGS would approach 50%.")
    print("   - Conclusion: This statement is NOT necessarily true.\n")

    # Statement D: The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.
    print("D. Epigenetic effects would limit h2 to less than 0.5.")
    print("   - Logic: Epigenetic effects are part of the environmental variance (Ve). The premise H2=0.5 already implies Ve=0.5.")
    print("   - The presence of environmental variance does not require the presence of non-additive genetic variance (Vd or Vi).")
    print(f"   - Counterexample: We can have Ve={Vp-Vg} (containing epigenetic effects) while Vg is purely additive (Va={Vg}), making h2 = H2 = {H2}. It is not strictly less than 0.5.")
    print("   - Conclusion: This statement is NOT necessarily true.\n")

analyze_heritability()