import sys

def analyze_heritability():
    """
    Analyzes statements about heritability and polygenic scores based on given parameters.
    """
    # --- Setup and Definitions ---
    # The problem states that for a phenotype, the broad-sense heritability (H2) is 0.5.
    # Broad-sense heritability H2 = Vg / Vp, where Vg is total genetic variance and Vp is total phenotypic variance.
    # Total genetic variance Vg = Va + Vd + Vi (Additive + Dominance + Epistatic variance).
    # A standard Polygenic Score (PGS) is built from additive effects, so its maximum theoretical
    # variance explained is the narrow-sense heritability (h2).
    # Narrow-sense heritability h2 = Va / Vp.

    H2 = 0.5
    print(f"Given: Broad-sense heritability (H2) = Vg / Vp = {H2}")
    print("A standard PGS can, at its best, explain h2 = Va / Vp of the variance.\n")

    # To demonstrate, let's assume a total phenotypic variance (Vp) of 100 units.
    Vp = 100
    # From H2 = 0.5, we can calculate the total genetic variance (Vg).
    Vg = H2 * Vp
    print(f"Let's assume a hypothetical Vp = {Vp}. Then the total genetic variance Vg = {H2} * {Vp} = {Vg}.\n")

    # --- Statement A Analysis ---
    # "The polygenic score can not explain more than 50% of the variance in the phenotype"
    print("--- Analysis of Statement A ---")
    print("The variance explained by a PGS is based on genetic effects.")
    print("The TOTAL variance from ALL genetic effects (Vg) is 50% of the phenotypic variance (Vp).")
    print("A PGS, which models some or all of these genetic effects, cannot explain more variance than what is available.")
    print("Mathematically, Va (additive variance) is a component of Vg (total genetic variance), so Va <= Vg.")
    print(f"Therefore, the PGS variance explained (h2 = Va/Vp) must be less than or equal to H2 (Vg/Vp).")
    print(f"PGS R-squared <= {H2}")
    print("Conclusion: Statement A is necessarily TRUE.\n")

    # --- Statement B & C Analysis ---
    # B: "Given an arbitrarily large GWAS, the PGS will approach a variance explained of 50%"
    # C: "...the PGS will not approach a variance explained of 50% due to non-linear effects"
    print("--- Analysis of Statements B and C ---")
    print("These statements depend on whether all genetic variance is additive or not.")

    # Scenario 1: All genetic variance is purely additive (Vd = 0, Vi = 0).
    print("\nScenario 1: Assume all genetic variance is additive.")
    Va_1 = Vg
    Vd_1 = 0
    Vi_1 = 0
    h2_1 = Va_1 / Vp
    print(f"In this case, Va = Vg = {Vg}.")
    print(f"The variance explained by a perfect PGS would be h2 = Va / Vp = {Va_1} / {Vp} = {h2_1}")
    print("In this scenario, the PGS *does* approach 50% explained variance.")
    print("This shows that Statement C is NOT necessarily true, as we have found a valid counterexample.")

    # Scenario 2: Some genetic variance is non-additive (dominance or epistasis exists).
    print("\nScenario 2: Assume non-additive genetic variance exists.")
    Va_2 = 30  # Assume Va is less than Vg
    Vd_2 = 15
    Vi_2 = 5
    # The components must sum to Vg: Va_2 + Vd_2 + Vi_2 = 30 + 15 + 5 = 50. This is a valid scenario.
    h2_2 = Va_2 / Vp
    print(f"In this case, Va = {Va_2}, which is less than Vg = {Vg}.")
    print(f"The variance explained by a perfect PGS would be h2 = Va / Vp = {Va_2} / {Vp} = {h2_2}")
    print("In this scenario, the PGS approaches 30%, which is NOT 50%.")
    print("This shows that Statement B is NOT necessarily true, as we have found a valid counterexample.")
    print("Since neither B nor C is true in all valid scenarios, neither is *necessarily* true.\n")

    # --- Statement D Analysis ---
    # "The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5"
    print("--- Analysis of Statement D ---")
    print("Epigenetic effects are typically considered part of the environmental variance (Ve) or a separate component of Vp.")
    print("The existence of epigenetic variance does NOT necessitate the existence of non-additive *genetic* variance (Vd or Vi).")
    print("Let's revisit Scenario 1, where all genetic variance is additive (Va = Vg).")
    h2_D = Va_1 / Vp
    print(f"In that scenario, h2 = Va / Vp = {Va_1} / {Vp} = {h2_D}")
    print("h2 can equal 0.5 even if epigenetic effects exist and contribute to the other 50% of phenotypic variance.")
    print("The statement claims h2 would be strictly *less than* 0.5, which is false in this case.")
    print("Conclusion: Statement D is NOT necessarily true.\n")

if __name__ == '__main__':
    analyze_heritability()