def analyze_heritability_problem():
    """
    Analyzes the quantitative genetics problem to determine the necessarily true statement.
    This function uses print statements to explain the logic based on fundamental principles.
    """
    H_squared = 0.5
    
    print("--- Problem Analysis ---")
    print(f"We are given a broad-sense heritability (H²) of {H_squared}.")
    print("\nLet's define the key concepts:")
    print("Vp = Total Phenotypic Variance")
    print("Vg = Total Genetic Variance")
    print("Ve = Environmental Variance")
    print("So, Vp = Vg + Ve")
    
    print("\nGenetic Variance (Vg) can be partitioned:")
    print("Vg = Va + Vd + Vi, where:")
    print("  Va = Additive Genetic Variance (sum of linear effects)")
    print("  Vd = Dominance Variance (non-linear interactions at one locus)")
    print("  Vi = Epistatic Variance (non-linear interactions between loci)")

    print("\nHeritability definitions:")
    print(f"Broad-sense (H²) = Vg / Vp. We are given H² = {H_squared}.")
    print("Narrow-sense (h²) = Va / Vp.")
    
    print("\nPolygenic Score (PGS) and Variance Explained:")
    print("A standard PGS is built from GWAS, which estimates additive effects.")
    print("Therefore, the maximum variance a PGS can explain (R²_PGS) is the narrow-sense heritability (h²).")
    print("R²_PGS <= h²")
    
    print("\n--- Core Relationship ---")
    print("Because Vg = Va + Vd + Vi, and Vd and Vi are non-negative, it must be true that: Va <= Vg.")
    print("Dividing by Vp gives: (Va / Vp) <= (Vg / Vp).")
    print(f"This means: h² <= H².")
    print(f"So, for our problem: h² <= {H_squared}.")
    
    print("\n--- Evaluating the Statements ---")
    
    print("\n[A] The PGS can not explain more than 50% of the variance.")
    print("This means: R²_PGS <= 0.5.")
    print(f"We know R²_PGS <= h² and h² <= H² = {H_squared}.")
    print(f"Therefore, it must be true that R²_PGS <= {H_squared}.")
    print("Statement A is NECESSARILY TRUE.")
    
    print("\n[B] The PGS will approach a variance explained of 50%.")
    print("This implies h² will be equal to 0.5.")
    print("This is only true if h² = H², which means all genetic variance is additive (Vd=0 and Vi=0).")
    print("This condition is not guaranteed by the problem statement.")
    print("Statement B is NOT necessarily true.")
    
    print("\n[C] The PGS will not approach a variance explained of 50% due to non-linear effects.")
    print("This implies that h² must be strictly less than 0.5, because non-additive effects (Vd or Vi) must exist.")
    print("However, the problem does not forbid the theoretical possibility that all genetic variance is purely additive (Vd=0 and Vi=0), in which case h² would equal 0.5.")
    print("Statement C is NOT necessarily true.")

    print("\n[D] The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("Epigenetic effects contribute to the total phenotypic variance Vp, usually as part of environmental variance Ve.")
    print(f"The ratio Vg/Vp is fixed at {H_squared}. Epigenetic effects do not necessarily alter the ratio of Va/Vg.")
    print(f"It's still possible that Vg is entirely additive (Va = Vg), which would make h² = H² = {H_squared}.")
    print("Statement D is NOT necessarily true.")

    print("\n--- Conclusion ---")
    print("Only statement A is a logical necessity based on the information provided.")


if __name__ == '__main__':
    analyze_heritability_problem()
    print("\n<<<A>>>")