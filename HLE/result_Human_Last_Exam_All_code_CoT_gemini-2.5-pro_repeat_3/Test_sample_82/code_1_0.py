def analyze_genetics_problem():
    """
    Analyzes the provided genetics problem statement by statement using established definitions.
    """

    # 1. Define knowns and fundamental relationships
    H_squared = 0.5  # Broad-sense heritability (Vg / Vp)

    print("--- Problem Analysis ---")
    print("Key Definitions:")
    print("Vp = Phenotypic Variance")
    print("Vg = Genetic Variance = Va (Additive) + Vd (Dominance) + Vi (Epistasis)")
    print("Ve = Environmental Variance")
    print("H^2 (Broad-sense Heritability) = Vg / Vp")
    print("h^2 (Narrow-sense Heritability) = Va / Vp")
    print("A standard Polygenic Score (PGS) explains at most h^2 of the variance.")
    print("\nGiven Information:")
    print(f"H^2 = Vg / Vp = {H_squared}")

    print("\nCore Deduction:")
    print("Because Vg = Va + Vd + Vi, and variance cannot be negative, it must be that Va <= Vg.")
    print("Dividing by Vp gives: (Va / Vp) <= (Vg / Vp)")
    print(f"This means h^2 <= H^2. Therefore, h^2 <= {H_squared}.")
    print("The maximum variance a PGS can explain is h^2, which must be less than or equal to 50%.")
    print("-" * 25)

    # 2. Evaluate each statement
    print("\nStatement Evaluations:")

    # Statement A
    print("\n[A] The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print(f"  - The max variance a PGS can explain is h^2.")
    print(f"  - We deduced that h^2 <= {H_squared}.")
    print(f"  - Therefore, the PGS variance explained must be <= {H_squared}. It cannot be *more* than 50%.")
    print("  - Verdict: NECESSARILY TRUE.")

    # Statement B
    print("\n[B] Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print(f"  - This implies that h^2 would approach {H_squared}.")
    print(f"  - This only happens if Vd and Vi are zero, so that Va = Vg.")
    print("  - The problem does not state that non-additive effects are zero.")
    print("  - Verdict: NOT NECESSARILY TRUE.")

    # Statement C
    print("\n[C] ...the polygenic score ... will not approach ... 50% due to gene-gene interactions and other non-linear effects...")
    print("  - This implies that non-additive effects (Vd, Vi) *must* be greater than zero.")
    print("  - While biologically likely, the problem does not forbid the theoretical case where Vd=0 and Vi=0.")
    print("  - If Vd=Vi=0, then h^2 = H^2, and this statement would be false.")
    print("  - Verdict: NOT NECESSARILY TRUE.")

    # Statement D
    print("\n[D] The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("  - Epigenetic effects are part of the environmental variance (Ve), not the genetic variance (Vg).")
    print(f"  - The relationship h^2 <= {H_squared} is independent of the composition of Ve.")
    print(f"  - It is possible that h^2 = {H_squared} (if Vd=Vi=0) even if epigenetic effects exist.")
    print("  - The statement requires h^2 to be strictly *less than* 0.5, which is not guaranteed.")
    print("  - Verdict: NOT NECESSARILY TRUE.")


# Execute the analysis
analyze_genetics_problem()
print("\n<<<A>>>")