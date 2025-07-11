def analyze_heritability_statements():
    """
    Analyzes statements about heritability and polygenic scores based on given information.
    """

    # --- 1. Problem Setup ---
    # We are given a broad-sense heritability (H^2) of 0.5.
    # Let's use a hypothetical total phenotypic variance (Vp) of 100 for easy calculation.
    H2 = 0.5
    Vp = 100.0
    
    # From H^2 = Vg / Vp, we can calculate the total genetic variance (Vg).
    Vg = H2 * Vp

    print("--- Step 1: Defining the Parameters ---")
    print("Broad-sense heritability, H² = Vg / Vp")
    print(f"Given H² = {H2} and assuming Vp = {Vp}, total genetic variance Vg is:")
    print(f"Vg = {H2} * {Vp} = {Vg}\n")

    # --- 2. Statement-by-Statement Analysis ---

    # Statement A Analysis
    print("--- Analysis of Statement A ---")
    print("Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print("Reasoning: A polygenic score is a predictor based on genetic information.")
    print("The maximum possible variance that can be explained by *any* form of genetic information (including non-linear effects) is the total genetic variance, Vg.")
    print(f"The proportion of phenotypic variance explained by Vg is H², which is {H2} (or 50%).")
    print("Therefore, no PGS can explain more than 50% of the variance.")
    print("Conclusion: Statement A is TRUE.\n")

    # Statements B and C Analysis
    print("--- Analysis of Statements B and C ---")
    print("Statement B: Given an arbitrarily large GWAS, the PGS will approach a variance explained of 50%.")
    print("Statement C: Given an arbitrarily large GWAS, the PGS ... will not approach ... 50% due to gene-gene interactions and other non-linear effects.")
    print("\nReasoning: Standard PGS are built by linearly summing effects from GWAS.")
    print("Such a score's predictive power is limited by narrow-sense heritability, h² = Va / Vp.")
    print("Total genetic variance Vg is composed of: Vg = Va (additive) + Vd (dominance) + Vi (epistasis).")
    print("So, Va <= Vg, which means h² <= H².")
    
    # Consider a biologically plausible scenario for a "polygenic trait" where non-additive effects exist.
    print("\nLet's model a realistic scenario for a complex 'polygenic trait':")
    # Assign plausible values where Vd + Vi > 0
    Va_scenario = 35.0
    Vd_scenario = 10.0
    Vi_scenario = 5.0
    # The sum must equal the total genetic variance, Vg
    # 35.0 + 10.0 + 5.0 = 50.0
    
    # Calculate the narrow-sense heritability for this scenario
    h2_scenario = Va_scenario / Vp
    
    print(f"Assume Vg ({Vg}) is composed of Va={Va_scenario}, Vd={Vd_scenario}, and Vi={Vi_scenario}.")
    print("The theoretical maximum variance explained by a linear PGS (h²) would be:")
    print("h² = Va / Vp")
    print(f"   = {Va_scenario} / {Vp} = {h2_scenario}")
    
    print(f"\nIn this realistic case, the PGS variance explained ({h2_scenario*100}%) is less than H² ({H2*100}%).")
    print("Because it is standard to assume polygenic traits have non-additive components (dominance, epistasis), a linear PGS will not capture all genetic variance.")
    print("Conclusion: Statement B is FALSE. Statement C is TRUE.\n")

    # Statement D Analysis
    print("--- Analysis of Statement D ---")
    print("Statement D: The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("Reasoning: We already know h² <= H² by definition. Since H² = 0.5, it follows that h² <= 0.5.")
    print("This inequality holds regardless of epigenetics.")
    print("h² is strictly less than 0.5 if non-additive genetic variance (Vd or Vi) is greater than zero, not because of epigenetics.")
    print("Epigenetic factors are generally modeled as environmental or gene-environment interactions, they do not inherently force Va < Vg.")
    print("Conclusion: Statement D is NOT necessarily true.\n")

    # Final Summary
    print("--- Final Conclusion ---")
    print("Statements necessarily true (under standard assumptions for polygenic traits): A and C.")

if __name__ == '__main__':
    analyze_heritability_statements()
    # The combination of A and C being correct corresponds to choice E.
    print("<<<E>>>")