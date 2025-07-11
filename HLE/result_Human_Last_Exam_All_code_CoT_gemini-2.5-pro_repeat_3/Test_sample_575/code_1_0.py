def analyze_heritability_statements():
    """
    Analyzes statements about heritability and polygenic scores based on given information.
    """

    # --- Step 1: Define Given Information and Principles ---

    # We are given the broad-sense heritability (H^2).
    H_sq = 0.5

    # Broad-sense heritability (H^2) is the proportion of total phenotypic variance (Vp)
    # that is due to all genetic variance (Vg).
    # The fundamental equation is: H^2 = Vg / Vp
    #
    # Vg itself is composed of:
    # Vg = Va + Vd + Vi
    # Va = Additive genetic variance
    # Vd = Dominance genetic variance
    # Vi = Epistatic (gene-gene interaction) variance
    #
    # A standard Polygenic Score (PGS) from GWAS primarily captures the additive effects,
    # and its predictive power relates to narrow-sense heritability (h^2 = Va / Vp).

    print("--- Analysis of the Statements ---")
    print(f"Given: Broad-sense heritability (H^2) is {H_sq}.")
    print("This means the total contribution of all genetic factors (Vg) to the phenotype's variance (Vp) is exactly 50%.")
    print("The final equation is: Vg / Vp = 0.5")
    print("-" * 40)

    # --- Step 2: Evaluate Each Statement ---

    # A. The polygenic score can not explain more than 50% of the variance in the phenotype.
    print("\n[Analysis of A]")
    print("A PGS is a predictor based entirely on genetic information.")
    print("Therefore, the absolute maximum variance it can possibly explain is the total genetic variance, Vg.")
    print(f"The maximum proportion of variance that can be explained is Vg / Vp, which is defined as H^2 = {H_sq}.")
    print("Therefore, no genetic-based score can explain more than 50% of the variance.")
    print("Conclusion: Statement A is NECESSARILY TRUE.")

    # B. Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.
    print("\n[Analysis of B]")
    print("A standard PGS captures additive variance (Va). Its potential variance explained is h^2 = Va / Vp.")
    print(f"We only know that H^2 = (Va + Vd + Vi) / Vp = {H_sq}.")
    print("If there is any non-additive genetic variance (Vd > 0 or Vi > 0), then Va < Vg, and h^2 < H^2.")
    print("Since we cannot assume Vd and Vi are zero, we cannot say h^2 will necessarily approach 0.5.")
    print("Conclusion: Statement B is NOT necessarily true.")

    # C. Given an arbitrarily large GWAS, the polygenic score ... will not approach a variance explained of 50%...
    print("\n[Analysis of C]")
    print("This statement assumes that non-additive effects (Vd or Vi) MUST exist.")
    print("However, it's theoretically possible for a trait that all genetic variance is additive (Vg = Va).")
    print(f"In that case, h^2 would equal H^2, and the variance explained would approach {H_sq}.")
    print("Since we can't rule out this possibility, the statement is not necessarily true.")
    print("Conclusion: Statement C is NOT necessarily true.")

    # D. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.
    print("\n[Analysis of D]")
    print("By definition, narrow-sense heritability (h^2) must be less than or equal to broad-sense heritability (H^2).")
    print(f"So, h^2 is already <= {H_sq}.")
    print("This statement claims h^2 must be strictly LESS than 0.5 if epigenetic effects exist.")
    print("Epigenetic effects are typically considered part of environmental variance (Ve), not genetic variance (Vg).")
    print("Their existence does not force Va to be less than Vg. It's possible that Vg = Va and h^2 = 0.5.")
    print("Conclusion: Statement D is NOT necessarily true.")

# Run the analysis
analyze_heritability_statements()

# --- Final Answer ---
# Based on the logical deduction, statement A is the only one that must be true.
print("\n<<<A>>>")