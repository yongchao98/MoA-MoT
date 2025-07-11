def analyze_heritability():
    """
    Analyzes the statements about heritability and polygenic scores based on the problem's premise.
    """
    # 1. Define parameters from the problem
    # H2 is broad-sense heritability. Vp is total phenotypic variance.
    H2 = 0.5
    Vp = 1.0  # Assume total variance is 1.0 for simplicity

    # From H2 = Vg / Vp, we can calculate total genetic variance (Vg)
    Vg = H2 * Vp
    # The remaining variance is environmental (Ve)
    Ve = Vp - Vg

    print("--- Initial Problem Setup ---")
    print(f"Broad-sense heritability (H^2) = {H2}")
    print(f"Assuming Total Phenotypic Variance (Vp) = {Vp}")
    print(f"Total Genetic Variance (Vg = H^2 * Vp) = {Vg}")
    print(f"Environmental Variance (Ve = Vp - Vg) = {Ve}\n")
    print("An ideal Polygenic Score (PGS) from a very large GWAS explains a proportion of variance equal to the narrow-sense heritability (h^2).\n")

    # 2. Define two possible scenarios for the composition of Genetic Variance (Vg)

    # Scenario 1: All genetic variance is purely additive.
    print("--- Scenario 1: Purely Additive Genetic Model ---")
    Va1 = Vg
    Vd1 = 0.0
    Vi1 = 0.0
    # Calculate narrow-sense heritability (h^2 = Va / Vp)
    h2_1 = Va1 / Vp
    R2_PGS_1 = h2_1
    print(f"Additive Variance (Va) = {Va1}")
    print(f"Non-Additive Variance (Vd + Vi) = {Vd1 + Vi1}")
    print(f"Narrow-sense heritability (h^2) = Va / Vp = {Va1} / {Vp} = {h2_1}")
    print(f"Ideal PGS Variance Explained = {R2_PGS_1}\n")

    # Scenario 2: Genetic variance is a mix of additive and non-additive effects.
    print("--- Scenario 2: Mixed Genetic Model (with non-additive effects) ---")
    Va2 = 0.6 * Vg  # Assume 60% of Vg is additive
    Vd2 = 0.2 * Vg  # 20% is dominance
    Vi2 = 0.2 * Vg  # 20% is epistasis
    h2_2 = Va2 / Vp
    R2_PGS_2 = h2_2
    print(f"Additive Variance (Va) = {Va2:.2f}")
    print(f"Non-Additive Variance (Vd + Vi) = {Vd2 + Vi2:.2f}")
    print(f"Narrow-sense heritability (h^2) = Va / Vp = {Va2:.2f} / {Vp} = {h2_2}")
    print(f"Ideal PGS Variance Explained = {R2_PGS_2}\n")

    # 3. Evaluate each statement using the scenarios
    print("--- Evaluating the Statements ---")

    # Statement A
    print("A. The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print(f"In Scenario 1, PGS explains {R2_PGS_1*100}%. In Scenario 2, it explains {R2_PGS_2*100}%.")
    print("In all possible cases, the PGS variance explained (h^2) must be less than or equal to the total genetic variance (H^2), which is 50%.")
    print("--> Statement A is NECESSARILY TRUE.\n")

    # Statement B
    print("B. Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print(f"This is only true in Scenario 1. It is false in Scenario 2, where the PGS only approaches {R2_PGS_2*100}%.")
    print("--> Statement B is NOT NECESSARILY TRUE.\n")

    # Statement C
    print("C. ...the PGS... will not approach... 50% due to gene-gene interactions and other non-linear effects...")
    print("This statement correctly identifies that a linear PGS only captures additive variance (Va).")
    print("If non-linear effects exist (like in Scenario 2), they contribute to H^2 but not h^2, so the PGS will fall short of H^2 (50%).")
    print("The statement describes a fundamental limitation of the tool. This limitation is always true.")
    print("--> Statement C is NECESSARILY TRUE.\n")

    # Statement D
    print("D. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print(f"Epigenetic effects contribute to Environmental Variance (Ve = {Ve}).")
    print(f"In Scenario 1, we have Ve > 0, but the narrow-sense heritability (h^2) is exactly {h2_1}.")
    print("The presence of non-genetic variance does not force the genetic variance to be non-additive.")
    print("--> Statement D is NOT NECESSARILY TRUE.\n")

    print("--- Final Conclusion ---")
    print("Statements A and C are the only ones that are necessarily true.")

analyze_heritability()