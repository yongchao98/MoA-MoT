def analyze_heritability():
    """
    Analyzes genetic concepts based on the problem statement to determine the correct answer.
    """
    # 1. Define initial parameters from the problem
    Vp = 1.0  # Let's define total phenotypic variance as 1.0 for simplicity
    H2 = 0.5  # Broad-sense heritability (given)

    # 2. Calculate the components of variance based on the premise
    Vg = H2 * Vp  # Total genetic variance
    Ve = Vp - Vg   # Environmental variance

    print("--- Foundational Concepts ---")
    print(f"Based on the premise (H2 = {H2}), we can set up the variance components:")
    print(f"Total Phenotypic Variance (Vp) = {Vp:.2f}")
    print(f"Broad-Sense Heritability (H2) = Vg / Vp = {Vg:.2f} / {Vp:.2f} = {H2:.2f}")
    print(f"This means the total Genetic Variance (Vg) is {Vg:.2f}, and Environmental Variance (Ve) is {Ve:.2f}.")
    print("-" * 35)

    # 3. Model two possible scenarios for the composition of Genetic Variance (Vg)
    print("\n--- Scenario 1: All Genetic Variance is Additive ---")
    Va_1 = Vg
    Vd_1 = 0.0
    Vi_1 = 0.0
    h2_1 = Va_1 / Vp
    print(f"Here, we assume no dominance or epistasis.")
    print(f"Vg = Va + Vd + Vi  ==>  {Vg:.2f} = {Va_1:.2f} + {Vd_1:.2f} + {Vi_1:.2f}")
    print(f"Narrow-Sense Heritability (h2) = Va / Vp = {Va_1:.2f} / {Vp:.2f} = {h2_1:.2f}")
    print(f"A perfect Polygenic Score (PGS) would explain {h2_1*100:.0f}% of variance.")

    print("\n--- Scenario 2: Genetic Variance includes Non-Additive Effects ---")
    # Assume additive variance is 60% of total genetic variance
    Va_2 = 0.6 * Vg
    Vd_2 = 0.2 * Vg
    Vi_2 = 0.2 * Vg
    h2_2 = Va_2 / Vp
    print(f"Here, we assume dominance and epistasis exist.")
    print(f"Vg = Va + Vd + Vi  ==>  {Vg:.2f} = {Va_2:.2f} + {Vd_2:.2f} + {Vi_2:.2f}")
    print(f"Narrow-Sense Heritability (h2) = Va / Vp = {Va_2:.2f} / {Vp:.2f} = {h2_2:.2f}")
    print(f"A perfect Polygenic Score (PGS) would explain {h2_2*100:.0f}% of variance.")
    print("-" * 35)

    # 4. Evaluate each statement
    print("\n--- Evaluating the Statements ---")
    # Statement A: PGS cannot explain more than 50% of variance.
    # The max variance a PGS can explain is H2 = 0.5.
    print(f"A: Can PGS explain > 50%? In Scenario 1, it explains {h2_1*100:.0f}%. In Scenario 2, {h2_2*100:.0f}%.")
    print("   Neither is > 50%. The max possible is H2 (50%). This statement is NECESSARILY TRUE.")

    # Statement B: PGS will approach 50% explained variance.
    print(f"B: Will PGS approach 50%? This is true for Scenario 1 ({h2_1*100:.0f}%) but false for Scenario 2 ({h2_2*100:.0f}%).")
    print("   Therefore, it is NOT necessarily true.")

    # Statement C: PGS will not approach 50% due to non-linear effects.
    print(f"C: Will PGS *not* approach 50%? This is true for Scenario 2 but false for Scenario 1.")
    print("   Therefore, it is NOT necessarily true.")
    
    # Statement D: Epigenetic effects would limit h2 to less than 0.5.
    print(f"D: Do epigenetic effects force h2 < 0.5? Epigenetic effects are part of Ve ({Ve:.2f}).")
    print(f"   In Scenario 1, even with Ve > 0, h2 equals {h2_1:.2f}. The statement is false.")
    print("   Therefore, it is NOT necessarily true.")

analyze_heritability()