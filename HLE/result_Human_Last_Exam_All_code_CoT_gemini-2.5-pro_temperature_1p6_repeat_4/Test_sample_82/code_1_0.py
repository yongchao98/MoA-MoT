def analyze_heritability_statements():
    """
    Analyzes the statements based on principles of quantitative genetics.
    """
    # Premise from the question
    H_squared = 0.5  # Broad-sense heritability (Vg / Vp)

    print("--- Problem Setup ---")
    print(f"Given: Broad-sense heritability (H^2) = Vg / Vp = {H_squared}")
    print("Definitions:")
    print("  Vp: Total Phenotypic Variance")
    print("  Vg: Total Genetic Variance (Vg = Va + Vd + Vi)")
    print("  Va: Additive Genetic Variance")
    print("  Vd: Dominance Variance (non-linear)")
    print("  Vi: Epistatic Variance (non-linear, gene-gene interaction)")
    print("  h^2: Narrow-sense heritability = Va / Vp")
    print("  A standard Polygenic Score (PGS) from GWAS explains variance approaching h^2.\n")

    # --- Analysis of Statement A ---
    print("--- Analysis of Statement A ---")
    print("Statement A: The PGS cannot explain more than 50% of the variance in the phenotype.")
    print("The maximum possible variance that ANY genetic predictor can explain is the total genetic variance, Vg.")
    print(f"The proportion of variance from all genetic factors is H^2 = Vg / Vp = {H_squared}.")
    print("The PGS explains h^2 = Va / Vp. Since Va is a component of Vg, Va <= Vg.")
    print(f"Therefore, h^2 <= H^2. This means the variance explained by a PGS is always less than or equal to {H_squared}.")
    print("Conclusion: Statement A is NECESSARILY TRUE.\n")

    # --- Analysis of Statements B and C ---
    print("--- Analysis of Statements B and C ---")
    print("Statement C claims the PGS will NOT approach 50% due to non-linear effects (dominance, epistasis).")
    print("This hinges on whether h^2 is strictly less than H^2.")
    print("\nLet's model a realistic scenario for a polygenic trait. Assume Vp = 100 units.")
    Vp = 100.0
    Vg = H_squared * Vp
    # Let's partition Vg into additive and non-additive components
    Va = 35.0
    Vd = 10.0
    Vi = 5.0
    print(f"Let Vp = {Vp}. Then total genetic variance Vg = {H_squared} * {Vp} = {Vg}.")
    print(f"Let's assume Vg is composed of: Va={Va}, Vd={Vd}, Vi={Vi}.")
    print(f"Check: {Va} + {Vd} + {Vi} = {Va + Vd + Vi}, which equals Vg.")
    h_squared = Va / Vp
    print(f"In this scenario, narrow-sense heritability h^2 = Va / Vp = {Va} / {Vp} = {h_squared}.")
    print(f"An ideal PGS would approach explaining {h_squared*100}% of the variance.")
    print(f"This is less than {H_squared*100}% precisely because the PGS misses the variance from Vd and Vi ({Vd+Vi} units).")
    print("This scenario aligns with the phrasing of the question, which distinguishes 'broad-sense' heritability.")
    print("Conclusion: Statement C is TRUE under this standard interpretation.\n")

    # --- Analysis of Statement D ---
    print("--- Analysis of Statement D ---")
    print("Statement D: The existence of any epigenetic effects would limit h^2 to less than 0.5.")
    print("Let's test this with a counterexample.")
    # Define a scenario where all genetic variance is additive
    Va_d = 50.0
    Vd_d = 0.0
    Vi_d = 0.0
    Vg_d = Va_d + Vd_d + Vi_d
    # Define some variance from epigenetics
    V_epi = 50.0
    # The new total phenotypic variance
    Vp_d = Vg_d + V_epi
    print(f"Assume a case where all genetic variance is additive: Vg = Va = {Vg_d}.")
    print(f"And assume epigenetic variance exists: V_epi = {V_epi}.")
    print(f"The total phenotypic variance is Vp = Vg + V_epi = {Vg_d} + {V_epi} = {Vp_d}.")
    H_squared_d = Vg_d / Vp_d
    h_squared_d = Va_d / Vp_d
    print(f"In this population, H^2 = Vg / Vp = {Vg_d} / {Vp_d} = {H_squared_d}, matching the premise.")
    print(f"And h^2 = Va / Vp = {Va_d} / {Vp_d} = {h_squared_d}.")
    print("Here, epigenetic effects exist, but h^2 is EQUAL to 0.5, not less than 0.5.")
    print("Conclusion: Statement D is NOT necessarily true.\n")

    # --- Final Conclusion ---
    print("--- Final Summary ---")
    print("Based on the analysis:")
    print(" - Statement A is necessarily true.")
    print(" - Statement C is true under the standard interpretation of genetic terminology.")
    print(" - Statements B and D are not necessarily true.")
    print("The correct option combines A and C.")

if __name__ == '__main__':
    analyze_heritability_statements()