def analyze_heritability():
    """
    Analyzes the relationship between heritability components and polygenic scores.
    """
    # --- Definitions and Given Information ---
    # Vp: Total Phenotypic Variance
    # Vg: Total Genetic Variance
    # Va: Additive Genetic Variance
    # Vd: Dominance Variance
    # Vi: Epistatic (Interaction) Variance
    #
    # H2 (Broad-sense heritability) = Vg / Vp
    # h2 (Narrow-sense heritability) = Va / Vp
    #
    # A standard Polygenic Score's (PGS) maximum explained variance is h2.
    
    H2 = 0.5  # Given broad-sense heritability

    # For demonstration, let's assume total phenotypic variance Vp = 1.0
    Vp = 1.0
    Vg = H2 * Vp

    print("--- Analysis of Heritability and Polygenic Scores ---")
    print(f"Given Broad-Sense Heritability (H2) = {H2}")
    print(f"This means Total Genetic Variance (Vg) is {H2*100}% of Phenotypic Variance (Vp).")
    print("\nThe maximum variance a standard PGS can explain is the Narrow-Sense Heritability (h2 = Va/Vp).")
    print("By definition, Vg = Va + Vd + Vi. Since Vd and Vi >= 0, it must be that Va <= Vg.")
    print("Therefore, h2 (Va/Vp) must be less than or equal to H2 (Vg/Vp).")
    print("-" * 60)

    # --- Illustrative Scenarios ---
    print("Let's consider a few scenarios where Vg = 0.5:")

    # Scenario 1: All genetic variance is additive
    Va1, Vd1, Vi1 = 0.5, 0.0, 0.0
    h2_1 = Va1 / Vp
    print("\nScenario 1: All genetic variance is additive (Vd=0, Vi=0)")
    print(f"  Va={Va1}, Vd={Vd1}, Vi={Vi1}  => Vg = {Va1+Vd1+Vi1}")
    print(f"  Max PGS Explained Variance (h2) = {h2_1:.2f} or {h2_1*100:.0f}%")
    print(f"  Is h2 <= H2?  {h2_1:.2f} <= {H2}  => {h2_1 <= H2}")

    # Scenario 2: Mix of additive and non-additive variance
    Va2, Vd2, Vi2 = 0.3, 0.1, 0.1
    h2_2 = Va2 / Vp
    print("\nScenario 2: Mix of additive and non-additive effects")
    print(f"  Va={Va2}, Vd={Vd2}, Vi={Vi2}  => Vg = {Va2+Vd2+Vi2}")
    print(f"  Max PGS Explained Variance (h2) = {h2_2:.2f} or {h2_2*100:.0f}%")
    print(f"  Is h2 <= H2?  {h2_2:.2f} <= {H2}  => {h2_2 <= H2}")

    # Scenario 3: Another mix
    Va3, Vd3, Vi3 = 0.2, 0.2, 0.1
    h2_3 = Va3 / Vp
    print("\nScenario 3: Additive variance is a smaller component")
    print(f"  Va={Va3}, Vd={Vd3}, Vi={Vi3}  => Vg = {Va3+Vd3+Vi3}")
    print(f"  Max PGS Explained Variance (h2) = {h2_3:.2f} or {h2_3*100:.0f}%")
    print(f"  Is h2 <= H2?  {h2_3:.2f} <= {H2}  => {h2_3 <= H2}")
    
    print("-" * 60)
    print("\nConclusion:")
    print("In all possible scenarios, the maximum variance a PGS can explain (h2) is less than or equal to the broad-sense heritability (H2=0.5).")
    print("Therefore, the polygenic score cannot explain MORE than 50% of the variance.")
    print("This makes statement A necessarily true.")

analyze_heritability()