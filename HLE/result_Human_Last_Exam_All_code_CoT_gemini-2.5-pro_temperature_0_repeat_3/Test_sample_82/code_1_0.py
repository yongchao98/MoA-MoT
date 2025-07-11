def analyze_heritability():
    """
    This function simulates population variance components to analyze the given statements.
    """
    # 1. Define baseline parameters from the problem
    Vp = 100.0  # Assume total phenotypic variance is 100 for clarity
    H2 = 0.5    # Given broad-sense heritability

    # 2. Calculate total genetic and environmental variance
    Vg = H2 * Vp
    Ve = Vp - Vg

    # 3. Decompose genetic variance into components for a typical polygenic trait
    # Vg = Va (additive) + Vd (dominance) + Vi (epistasis)
    # Assume non-additive effects exist, so Va < Vg
    Va = 35.0
    Vd = 10.0
    Vi = 5.0

    # 4. Calculate narrow-sense heritability (h^2)
    h2 = Va / Vp

    # 5. The max variance explained by a linear PGS is h^2
    max_r2_pgs = h2
    
    print("--- Analysis of Heritability and Polygenic Scores ---")
    print(f"Given H^2 = {H2}, we can model the variance components:")
    print(f"Total Phenotypic Variance (Vp) = {Vp}")
    print(f"Total Genetic Variance (Vg) = {Vg}")
    print(f"  - Additive Variance (Va)   = {Va}")
    print(f"  - Dominance Variance (Vd)  = {Vd}")
    print(f"  - Epistatic Variance (Vi)  = {Vi}")
    print(f"Environmental Variance (Ve)  = {Ve}")
    print("-" * 50)
    
    print("--- Evaluating the Statements ---")
    # Statement A: PGS cannot explain more than 50% of variance.
    # The ultimate ceiling for any genetic predictor is H^2.
    print(f"Statement A: The maximum variance any genetic predictor can explain is H^2 = {H2:.2f} (50%).")
    print("             Therefore, a PGS cannot explain more than 50%. This is TRUE.")
    
    # Statement C: Linear PGS will not approach 50% due to non-linear effects.
    # A linear PGS is limited by h^2. The gap between h^2 and H^2 is due to Vd and Vi.
    gap = H2 - h2
    gap_due_to_Vd_Vi = (Vd + Vi) / Vp
    print(f"\nStatement C: A linear PGS is limited by h^2 = Va/Vp = {Va}/{Vp} = {h2:.2f} (35%).")
    print(f"             It will not approach H^2 ({H2:.2f}) because it cannot capture non-linear effects.")
    print(f"             The variance missed is (Vd+Vi)/Vp = ({Vd}+{Vi})/{Vp} = {gap_due_to_Vd_Vi:.2f}. This is TRUE.")
    print("-" * 50)
    
    print("\nConclusion: Both statements A and C are necessarily true.")

analyze_heritability()