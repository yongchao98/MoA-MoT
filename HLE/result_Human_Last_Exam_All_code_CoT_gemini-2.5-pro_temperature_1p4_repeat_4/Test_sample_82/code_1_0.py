import sys

def analyze_heritability_statements():
    """
    Analyzes statements about heritability and polygenic scores based on given information.
    """
    # --- Given Information ---
    # We are told that the broad-sense heritability (H^2) is 0.5.
    H_squared = 0.5

    # --- Principles of Quantitative Genetics ---
    # Vp: Total Phenotypic Variance
    # Vg: Total Genetic Variance
    # Ve: Environmental Variance
    # Va: Additive Genetic Variance
    # Vd: Dominance Genetic Variance
    # Vi: Epistatic (Interaction) Genetic Variance
    #
    # Key relationships:
    # 1. Vp = Vg + Ve
    # 2. Vg = Va + Vd + Vi
    # 3. Broad-sense heritability H^2 = Vg / Vp
    # 4. Narrow-sense heritability h^2 = Va / Vp
    # 5. A Polygenic Score's (PGS) max explained variance is h^2.

    print("--- Problem Analysis ---")
    print(f"Given: Broad-sense heritability (H^2) = Vg / Vp = {H_squared}")
    print("A Polygenic Score (PGS) from GWAS explains variance due to additive effects, with its theoretical maximum being narrow-sense heritability (h^2).")
    print("\n")


    # --- Statement A Analysis ---
    print("--- Evaluating Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype. ---")
    print("1. From the definitions, Vg = Va + Vd + Vi.")
    print("2. Since variance components (Vd, Vi) cannot be negative, Va must be less than or equal to Vg (Va <= Vg).")
    print("3. Dividing by Vp gives us the relationship between the heritabilities: (Va / Vp) <= (Vg / Vp), which means h^2 <= H^2.")
    H_squared_val = 0.5
    print(f"4. With H^2 = {H_squared_val}, we have the inequality: h^2 <= {H_squared_val}.")
    print(f"5. Since the PGS can explain at most h^2 of the variance, it can not explain more than {H_squared_val*100}% of the variance.")
    print("Conclusion: Statement A is necessarily TRUE.\n")

    # --- Statement C Analysis ---
    print("--- Evaluating Statement C: Given an arbitrarily large GWAS, the PGS will not approach a variance explained of 50% due to non-linear effects... ---")
    print("1. This statement claims the variance explained by the PGS will be strictly less than 50%, i.e., h^2 < 0.5.")
    print("2. This is true only if the non-additive genetic variance (Vd + Vi) is greater than 0.")
    print("3. The terms 'polygenic score' and 'typical GWAS data' refer to complex traits. For any complex trait, it's a standard biological assumption that non-additive genetic effects exist (i.e., Vd + Vi > 0). The distinction between H^2 and h^2 is meaningful only in this context.")
    print("4. Under this standard assumption, Va is strictly less than Vg, which means h^2 is strictly less than H^2.")
    # We can show this with an example equation.
    # Let Vp = 100. Then from H^2 = 0.5, Vg = 50.
    # Assume non-additive effects exist. Let Va = 40, and (Vd + Vi) = 10.
    Va_example = 40
    non_additive_V_example = 10
    Vg_example = Va_example + non_additive_V_example
    h_squared_example_numerator = Va_example
    H_squared_example_numerator = Vg_example
    Vp_denominator = 100
    print(f"5. Example: If Va={h_squared_example_numerator} and (Vd+Vi)={non_additive_V_example}, then Vg = {Va_example} + {non_additive_V_example} = {Vg_example}.")
    print(f"   In this case, h^2 = {h_squared_example_numerator}/{Vp_denominator} = {h_squared_example_numerator/Vp_denominator}, which is less than H^2 = {H_squared_example_numerator}/{Vp_denominator} = {H_squared_example_numerator/Vp_denominator}.")
    print("6. Therefore, the PGS variance (approaching h^2) will not reach 50% because of these non-additive effects.")
    print("Conclusion: Statement C is TRUE under the standard interpretation of a polygenic trait.\n")

    # --- Final Summary ---
    print("--- Final Conclusion ---")
    print("Statement A is true from first principles.")
    print("Statement C is true based on the standard assumptions for a complex polygenic trait.")
    print("Statements B is the opposite of C, so it is false.")
    print("Statement D incorrectly links epigenetic effects to the limit on h^2.")
    print("Thus, both A and C are the correct statements.")


if __name__ == '__main__':
    analyze_heritability_statements()
<<<E>>>