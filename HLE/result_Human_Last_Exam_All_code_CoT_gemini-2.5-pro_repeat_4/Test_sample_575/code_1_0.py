def explain_heritability():
    """
    Demonstrates the relationship between heritability and Polygenic Scores (PGS).
    """
    # --- Setup as per the problem ---
    H2 = 0.5  # Broad-sense heritability

    # For demonstration, let's assume a total phenotypic variance (Vp)
    Vp = 100.0

    # --- Calculations based on definitions ---
    # Total genetic variance (Vg) is H2 * Vp
    Vg = H2 * Vp

    # --- Explanation ---
    print("--- Understanding Heritability and Polygenic Scores ---")
    print(f"Given Broad-Sense Heritability (H^2) = {H2}")
    print(f"This means that out of a total phenotypic variance (Vp) of {Vp},")
    print(f"the total genetic variance (Vg) is: Vg = {H2} * {Vp} = {Vg}\n")

    print("Total genetic variance (Vg) is composed of:")
    print("Vg = Va (Additive) + Vd (Dominance) + Vi (Epistasis)\n")

    print("A standard Polygenic Score (PGS) from GWAS primarily captures additive effects.")
    print("Its predictive power (R^2) is limited by the narrow-sense heritability (h^2).")
    print("h^2 = Va / Vp\n")

    print("The crucial relationship is that Additive Variance (Va) cannot be greater than Total Genetic Variance (Vg).")
    print("This means Va <= Vg.")
    print("Therefore, h^2 (which is Va/Vp) must be less than or equal to H^2 (which is Vg/Vp).\n")

    print("--- Let's test two possible scenarios ---\n")

    # --- Scenario 1: All genetic variance is additive ---
    print("Scenario 1: All genetic variance is additive (an extreme, but possible case).")
    Va_1 = Vg
    Vd_1 = 0.0
    Vi_1 = 0.0
    print(f"Equation: Vg = {Va_1} (Va) + {Vd_1} (Vd) + {Vi_1} (Vi) = {Va_1 + Vd_1 + Vi_1}")
    h2_1 = Va_1 / Vp
    print(f"Equation: Narrow-sense heritability (h^2) = Va / Vp = {Va_1} / {Vp} = {h2_1}")
    print(f"Maximum PGS R^2 in this scenario = {h2_1:.2%}\n")

    # --- Scenario 2: A mix of genetic variance components ---
    print("Scenario 2: Genetic variance includes non-additive effects (a more realistic case).")
    # Assign some values, ensuring they sum to Vg
    Va_2 = Vg * 0.7
    Vd_2 = Vg * 0.2
    Vi_2 = Vg * 0.1
    print(f"Equation: Vg = {Va_2:.1f} (Va) + {Vd_2:.1f} (Vd) + {Vi_2:.1f} (Vi) = {Va_2 + Vd_2 + Vi_2}")
    h2_2 = Va_2 / Vp
    print(f"Equation: Narrow-sense heritability (h^2) = Va / Vp = {Va_2:.1f} / {Vp} = {h2_2}")
    print(f"Maximum PGS R^2 in this scenario = {h2_2:.2%}\n")


    # --- Conclusion ---
    print("--- Conclusion ---")
    print(f"In Scenario 1, the max PGS R^2 is {h2_1:.0%}.")
    print(f"In Scenario 2, the max PGS R^2 is {h2_2:.0%}.")
    print("\nIn all possible cases, the variance explained by the PGS (max R^2 = h^2) is less than or equal to the broad-sense heritability (H^2 = 0.5).")
    print("The PGS can never explain MORE than 50% of the variance.")
    print("\nThis confirms that statement A is the only one that is necessarily true.")

if __name__ == "__main__":
    explain_heritability()