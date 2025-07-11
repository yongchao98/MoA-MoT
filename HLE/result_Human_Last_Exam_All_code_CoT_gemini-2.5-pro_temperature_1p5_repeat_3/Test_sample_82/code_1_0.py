def analyze_heritability_problem():
    """
    Analyzes the population genetics problem step-by-step and prints the reasoning.
    """

    # --- Setup from the problem statement ---
    H2 = 0.5  # Broad-sense heritability (Vg / Vp)

    print("--- Problem Analysis ---")
    print("We are given the following information:")
    print(f"Broad-sense Heritability (H^2) = Vg / Vp = {H2}")
    print("Where: Vp = Total Phenotypic Variance")
    print("       Vg = Total Genetic Variance")
    print("Genetic variance (Vg) can be decomposed into:")
    print("       Vg = Va + Vd + Vi")
    print("       Va = Additive Genetic Variance")
    print("       Vd = Dominance Variance")
    print("       Vi = Epistatic (Interaction) Variance")
    print("\nFrom this, we define Narrow-sense Heritability (h^2):")
    print("h^2 = Va / Vp")
    print("A standard Polygenic Score (PGS) from GWAS is an additive model,")
    print("so the maximum variance it can explain is h^2.")
    print("It is always true that Va <= Vg, therefore h^2 <= H^2.")
    print("-" * 25)
    print("\n--- Evaluating the Statements ---\n")

    # --- Statement A ---
    print("A. The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print("Reasoning:")
    print("The PGS is a predictor based on an individual's genetics.")
    print("The total variance explained by all genetic factors is Vg.")
    print(f"The maximum proportion of variance any genetic predictor can explain is Vg / Vp = H^2 = {H2}.")
    print(f"The variance explained by a PGS must be <= H^2.")
    print(f"Therefore, PGS explained variance <= {H2}.")
    print("Conclusion: Statement A is necessarily TRUE.\n")

    # --- Statement B ---
    print("B. Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print("Reasoning:")
    print(f"A standard PGS will approach explaining h^2 (Va / Vp), not H^2 ({H2}).")
    print(f"This statement would only be true if h^2 = H^2, which requires non-additive variance (Vd + Vi) to be 0.")
    print("Since the existence of non-additive effects is not ruled out, we cannot assume h^2 = H^2.")
    print("Conclusion: Statement B is NOT necessarily true.\n")
    
    # --- Statement C ---
    print("C. Given an arbitrarily large GWAS, the polygenic score constructed via linearly summing GWAS effect sizes")
    print("   will not approach a variance explained of 50% due to gene-gene interactions and other non-linear")
    print("   effects such as dominance.")
    print("Reasoning:")
    print("This statement claims the explained variance (h^2) will be less than H^2 (0.5).")
    print(f"This is true if there is any non-additive variance (Vd + Vi > 0), which would mean Va < Vg, and thus h^2 < H^2.")
    print("The term 'polygenic trait' implies a complex trait influenced by many genes, where interactions (epistasis) and dominance are the biological norm.")
    print("Under this standard interpretation, non-additive effects exist, and the additive PGS cannot capture the full genetic variance.")
    print(f"The shortfall between what the PGS captures (h^2) and the total genetic variance ({H2}) is precisely due to these non-linear effects.")
    print("Conclusion: Statement C is necessarily TRUE.\n")

    # --- Statement D ---
    print("D. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("Reasoning:")
    print(f"We already know that h^2 <= H^2, so h^2 is already limited to be less than or equal to {H2}.")
    print("This statement claims that epigenetic effects would force h^2 to be strictly less than 0.5.")
    print("Epigenetic effects are generally considered part of environmental variance (Ve), not genetic variance (Vg).")
    print("Their existence does not require the existence of non-additive GENETIC variance (Vd or Vi).")
    print(f"It is theoretically possible for h^2 = H^2 = {H2} even if epigenetic effects also exist.")
    print("Conclusion: Statement D is NOT necessarily true.\n")

    # --- Final Conclusion ---
    print("-" * 25)
    print("Summary:")
    print("- Statement A is TRUE.")
    print("- Statement C is TRUE.")
    print("Therefore, the correct choice combines A and C.")

if __name__ == '__main__':
    analyze_heritability_problem()
<<<E>>>