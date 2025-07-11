def analyze_heritability_statement():
    """
    Analyzes statements about heritability and polygenic scores based on given information.
    """

    # --- Step 1 & 2: Define concepts and their relationships ---

    # H_squared: Broad-sense heritability, the proportion of phenotypic variance (Vp)
    # due to all genetic variance (Vg).
    # H_squared = Vg / Vp
    H_squared = 0.5

    # Vg: Total genetic variance. It's composed of:
    # Va: Additive genetic variance (effects that sum up linearly).
    # V_non_additive: Non-additive variance (dominance and epistatic interactions).
    # Equation: Vg = Va + V_non_additive

    # R_squared_PGS: The variance in the phenotype explained by a Polygenic Score.
    # A standard PGS is based on summing linear effects from a GWAS.
    # Its theoretical maximum is the narrow-sense heritability (h_squared).
    # h_squared = Va / Vp
    # So, R_squared_PGS <= h_squared

    # --- Step 3: Analyze the statements ---

    print("Analyzing the problem based on genetic principles.")
    print(f"Given: Broad-sense heritability (H^2) = {H_squared}\n")

    print("The key relationships are:")
    print("1. Variance Explained by PGS (R_squared_PGS) <= Narrow-sense heritability (h^2)")
    print("2. Narrow-sense heritability (h^2) = Additive Genetic Variance (Va) / Phenotypic Variance (Vp)")
    print("3. Broad-sense heritability (H^2) = Total Genetic Variance (Vg) / Phenotypic Variance (Vp)")
    print("4. Total Genetic Variance (Vg) >= Additive Genetic Variance (Va)")
    print("From (2), (3), and (4), it follows that h^2 <= H^2.")
    print("Therefore, R_squared_PGS <= h^2 <= H^2.\n")


    print("--- Evaluating Each Statement ---\n")

    # Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype.
    print("Analysis of Statement A:")
    print(f"The analysis above shows that R_squared_PGS <= H^2.")
    print(f"Substituting the given value: R_squared_PGS <= {H_squared}.")
    print("This means the variance explained by a PGS cannot exceed 50%.")
    print("Conclusion: Statement A is NECESSARILY TRUE.\n")

    # Statement B: Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.
    print("Analysis of Statement B:")
    print("A perfect PGS (from an infinitely large GWAS) would explain h^2 of the variance.")
    print(f"This statement claims h^2 will approach {H_squared}.")
    print(f"This is only true if h^2 = H^2, which requires all genetic variance to be additive (V_non_additive = 0).")
    print("The problem does not guarantee this. Non-additive effects could exist, making h^2 < H^2.")
    print("Conclusion: Statement B is NOT necessarily true.\n")

    # Statement C: Given an arbitrarily large GWAS, the PGS ... will not approach ... 50% due to non-linear effects...
    print("Analysis of Statement C:")
    print("This statement claims that h^2 < H^2 is a necessity.")
    print("However, it is theoretically possible for a trait to have no non-additive genetic variance.")
    print(f"In such a case, V_non_additive = 0, meaning h^2 = H^2 = {H_squared}.")
    print("Since it's possible for h^2 to equal 0.5, we cannot state it will *not* approach 0.5.")
    print("Conclusion: Statement C is NOT necessarily true.\n")
    
    # Statement D: The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.
    print("Analysis of Statement D:")
    print("Epigenetic effects contribute to environmental variance (Ve), not genetic variance (Vg).")
    print(f"We are given H^2 = Vg / (Vg + Ve) = {H_squared}. This already accounts for Ve.")
    print(f"h^2 = Va / (Vg + Ve). The existence of epigenetic effects does not change the relationship between Va and Vg.")
    print("It is still possible for all genetic variance to be additive (Va = Vg), which would make h^2 = H^2 = 0.5.")
    print("Conclusion: Statement D is NOT necessarily true.\n")


if __name__ == '__main__':
    analyze_heritability_statement()
    print("<<<A>>>")
