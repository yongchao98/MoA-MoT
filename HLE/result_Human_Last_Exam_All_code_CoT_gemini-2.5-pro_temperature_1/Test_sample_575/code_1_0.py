def analyze_heritability_question():
    """
    Analyzes the relationship between heritability and polygenic scores
    to determine the correct statement among the given choices.
    """
    # --- Step 1: Define the given parameters and core concepts ---

    # For simplicity, we can normalize the total Phenotypic Variance (Vp) to 1.0.
    Vp = 1.0

    # The problem states that the broad-sense heritability (H_squared) is 0.5.
    H_squared = 0.5
    H_squared_percent = H_squared * 100

    # Broad-sense heritability is the proportion of phenotypic variance explained by
    # the total genetic variance (Vg).
    # Equation: H^2 = Vg / Vp
    Vg = H_squared * Vp

    # Total genetic variance (Vg) is the sum of additive (Va), dominance (Vd),
    # and epistatic (Vi) variance.
    # Equation: Vg = Va + Vd + Vi
    # From the problem, we know: Va + Vd + Vi = 0.5

    # Narrow-sense heritability (h_squared) is the proportion of phenotypic variance
    # explained by *only* the additive genetic variance (Va).
    # Equation: h^2 = Va / Vp

    # --- Step 2: Analyze the Polygenic Score (PGS) and its limits ---

    # A standard Polygenic Score (PGS) is built by linearly summing effects from a GWAS.
    # It primarily captures additive genetic effects. Therefore, its theoretical maximum
    # variance explained is equal to the narrow-sense heritability (h^2).
    # Max R_squared_PGS = h^2 = Va / Vp

    # Any possible predictor based on genetics (even one that could capture dominance
    # and epistasis perfectly) cannot explain more variance than the total genetic
    # variance (Vg).
    # Max R_squared_any_genetic_model = H^2 = Vg / Vp

    print("--- Problem Analysis ---")
    print(f"Broad-sense heritability (H^2) is given as {H_squared}.")
    print(f"This means the total genetic variance (Vg) accounts for {H_squared_percent:.0f}% of the total phenotypic variance (Vp).")
    print(f"Final Equation: Vg = H^2 * Vp = {H_squared} * {Vp} = {Vg}")
    print("\n--- Evaluating Answer Choices ---")

    # --- Step 3: Evaluate each statement ---

    # Choice A
    print("\n[A] The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print(f"The PGS is a genetic predictor. The absolute maximum variance that can be explained by *any* genetic factor is the broad-sense heritability, H^2, which is {H_squared_percent:.0f}%.")
    print("Therefore, no PGS, regardless of its construction, can explain more variance than this limit.")
    print("This statement is NECESSARILY TRUE.")

    # Choice B
    print("\n[B] Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print("A standard PGS can, at best, capture the additive variance (Va), explaining h^2 = Va / Vp of phenotypic variance.")
    print("This statement assumes that h^2 approaches H^2 (i.e., Va approaches Vg).")
    print("This would only happen if dominance (Vd) and epistatic (Vi) effects were zero. The problem does not state this.")
    print("If Vd > 0 or Vi > 0, then h^2 will be strictly less than H^2 (0.5).")
    print("This statement is NOT NECESSARILY TRUE.")

    # Choice C
    print("\n[C] Given an arbitrarily large GWAS, the PGS ... will not approach a variance explained of 50% due to ... non-linear effects.")
    print("This statement claims the opposite of B, that it's impossible for h^2 to equal H^2.")
    print("However, it is theoretically possible that all genetic variance is purely additive (Vd = 0 and Vi = 0).")
    print("In that specific, but possible, scenario, h^2 would equal H^2 (0.5), and the PGS would approach explaining 50% of the variance.")
    print("This statement is NOT NECESSARILY TRUE.")

    # Choice D
    print("\n[D] The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("Heritability partitions variance into genetic (Vg) and environmental (Ve) components. Epigenetic effects are typically modeled as part of Ve.")
    print(f"Here, Vg = {Vg} and Ve = {1.0 - Vg}. The existence of epigenetic effects (part of Ve) does not prevent h^2 from reaching its maximum possible value of 0.5 (which would happen if Vg = Va).")
    print("This statement is NOT NECESSARILY TRUE.")


if __name__ == '__main__':
    analyze_heritability_question()
