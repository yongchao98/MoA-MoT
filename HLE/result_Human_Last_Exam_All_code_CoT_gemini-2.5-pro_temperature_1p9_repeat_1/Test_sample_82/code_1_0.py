def analyze_heritability_statements():
    """
    Analyzes statements about heritability and polygenic scores
    based on the given parameters.
    """
    # --- Step 1: Set up parameters from the problem ---
    # Broad-sense heritability (H_squared) is given as 0.5.
    H_squared = 0.5
    # For demonstration, let's assume a total phenotypic variance (Vp) of 100 units.
    Vp = 100.0

    # --- Step 2: Calculate the main variance components based on H_squared ---
    # Vg is the total variance from all genetic sources.
    Vg = H_squared * Vp
    # Ve is the variance from environmental sources.
    Ve = Vp - Vg

    # --- Step 3: Model the sub-components of genetic variance (Vg) ---
    # For a typical complex (polygenic) trait, Vg = Va + Vd + Vi.
    # Non-additive effects (Vd, Vi) are generally present. We'll model this by
    # assuming that Va is a fraction of Vg, and Vd and Vi make up the rest.
    # A standard assumption for many traits is that Va is the largest component.
    Va = 0.35 * Vp  # Additive variance (35% of total variance)
    Vd = 0.10 * Vp  # Dominance variance (10% of total variance)
    Vi = 0.05 * Vp  # Epistatic variance (5% of total variance)
    # Check if the components sum to the total genetic variance
    # Vg = 35.0 + 10.0 + 5.0 = 50.0, which matches H_squared * Vp.

    # --- Step 4: Calculate narrow-sense heritability (h_squared) ---
    # h_squared is the proportion of phenotypic variance explained by ADDITIVE effects.
    # This is the theoretical maximum variance explainable by a standard linear PGS.
    h_squared = Va / Vp

    print("--- Quantitative Genetics Model ---")
    print(f"Given Broad-Sense Heritability (H²): {H_squared}")
    print(f"Assumed Phenotypic Variance (Vp): {Vp:.1f}")
    print("\nCalculated Variance Components:")
    print(f"Total Genetic Variance (Vg) = H² * Vp => {H_squared} * {Vp:.1f} = {Vg:.1f}")
    print(f"Environmental Variance (Ve) = Vp - Vg => {Vp:.1f} - {Vg:.1f} = {Ve:.1f}")
    print("\nAssumed breakdown of Genetic Variance (Vg):")
    print(f"Vg = Va + Vd + Vi => {Vg:.1f} = {Va:.1f} + {Vd:.1f} + {Vi:.1f}")
    print("\nCalculated Narrow-Sense Heritability (h²):")
    print(f"h² = Va / Vp => {Va:.1f} / {Vp:.1f} = {h_squared}")
    print("-------------------------------------------\n")

    print("--- Evaluating the Statements ---\n")

    # Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype.
    print("A: The polygenic score is based on genetics, so its explanatory power is capped by the total genetic variance (Vg).")
    print(f"The total genetic variance as a fraction of total variance is H², which is {H_squared:.0%}.")
    print("Therefore, no PGS, regardless of its design, can explain more than 50% of the variance.")
    print("Verdict: Statement A is NECESSARILY TRUE.\n")

    # Statement B: Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.
    print("B: A standard linear PGS captures additive genetic variance (Va), not total genetic variance (Vg).")
    print(f"Its theoretical maximum is h², which we calculated as {h_squared:.0%}.")
    print(f"It would only approach 50% if h² equaled H², which requires dominance (Vd) and epistatic (Vi) variance to be zero.")
    print("Verdict: Statement B is NOT necessarily true.\n")

    # Statement C: Given an arbitrarily large GWAS, the polygenic score constructed via linearly summing GWAS effect sizes
    # will not approach a variance explained of 50% due to gene-gene interactions and other non-linear effects such as dominance.
    print("C: This statement correctly identifies why the PGS performance is limited.")
    print(f"The total genetic variance H² is the sum of additive and non-additive components: H² = h² + (Vd + Vi) / Vp.")
    print(f"The equation for our model is: {H_squared} = {h_squared} + ({Vd:.1f} + {Vi:.1f}) / {Vp:.1f}")
    print("Because non-additive effects (Vd and Vi) exist for complex traits, h² is less than H².")
    print("A linear PGS approaches h², so it will not approach H² (50%).")
    print("Verdict: Statement C is NECESSARILY TRUE for any non-trivial polygenic trait.\n")

    # Statement D: The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.
    print("D: By definition, h² is less than or equal to H² (0.5). This is due to non-additive GENETIC variance (Vd, Vi).")
    print("Epigenetic effects that are not encoded in the DNA sequence are typically considered part of the environmental variance (Ve).")
    print("This statement misattributes the reason why h² < H².")
    print("Verdict: Statement D is NOT necessarily true.\n")

    print("-------------------------------------------")
    print("Conclusion: Both statements A and C are necessarily true.")
    print("The correct answer choice is E.")

if __name__ == '__main__':
    analyze_heritability_statements()