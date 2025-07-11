def analyze_heritability():
    """
    Illustrates the concepts of heritability and polygenic scores with a numerical example.
    """
    # Let's assume the total phenotypic variance (Vp) is 100 units for simplicity.
    Vp = 100.0

    # From the problem, broad-sense heritability (H^2) is 0.5.
    H_squared = 0.5

    # This means the total genetic variance (Vg) is H^2 * Vp.
    Vg = H_squared * Vp

    print("--- Setup based on the problem statement ---")
    print(f"Total Phenotypic Variance (Vp): {Vp}")
    print(f"Broad-Sense Heritability (H^2): {H_squared}")
    print(f"Total Genetic Variance (Vg = H^2 * Vp): {Vg}")
    print("-" * 40)

    # --- Analysis of Statement A ---
    print("--- Analysis of Statement A ---")
    print("A polygenic score (PGS) is a model of the genetic contribution to a phenotype.")
    print("Its predictive power cannot exceed the total variance caused by all genetic factors (Vg).")
    print(f"The maximum proportion of variance any genetic model can explain is Vg / Vp = {Vg} / {Vp} = {H_squared}.")
    print("Therefore, the PGS cannot explain more than 50% of the variance. Statement A is TRUE.")
    print("-" * 40)

    # --- Analysis of Statement C ---
    print("--- Analysis of Statement C ---")
    print("A standard linear PGS captures only additive genetic variance (Va).")
    print("Let's model a realistic polygenic trait with non-additive effects (dominance Vd, epistasis Vi).")
    
    # We set Va to be less than Vg, and allocate the rest to Vd and Vi.
    Va = 40.0
    Vd = 5.0
    Vi = 5.0
    
    # The sum must equal our total genetic variance Vg.
    Vg_check = Va + Vd + Vi
    
    print(f"Let's assume a breakdown of Vg ({Vg}): Va={Va}, Vd={Vd}, Vi={Vi}. (Check: {Va}+{Vd}+{Vi}={Vg_check})")

    # The variance explained by a perfect linear PGS approaches narrow-sense heritability (h^2).
    h_squared = Va / Vp

    print(f"\nThe narrow-sense heritability (h^2) in this scenario is Va / Vp.")
    # Here we output the numbers in the final equation as requested
    print(f"Final Equation: h^2 = {Va} / {Vp} = {h_squared}")
    print(f"The PGS approaches explaining {h_squared*100}% of the variance.")
    print(f"Since non-additive effects exist (Vd+Vi > 0), h^2 ({h_squared}) is less than H^2 ({H_squared}).")
    print("Thus, the linear PGS will not approach explaining 50% of the variance. Statement C is TRUE.")
    print("-" * 40)

# Run the analysis
analyze_heritability()