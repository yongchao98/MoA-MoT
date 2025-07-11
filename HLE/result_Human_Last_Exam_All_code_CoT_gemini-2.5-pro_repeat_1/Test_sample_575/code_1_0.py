import sys

def analyze_heritability_statements():
    """
    Analyzes statements about heritability and polygenic scores based on given parameters.
    """
    # --- Setup from the problem statement ---
    # We are given the broad-sense heritability (H2).
    H2 = 0.5

    # For demonstration purposes, let's assume a total phenotypic variance (Vp).
    # The specific value doesn't matter, but it makes the numbers concrete.
    Vp = 100.0

    # From the definition H2 = Vg / Vp, we can calculate the total genetic variance (Vg).
    Vg = H2 * Vp

    print("--- Problem Setup ---")
    print(f"Broad-sense heritability (H²): {H2}")
    print(f"Assumed total phenotypic variance (Vp): {Vp}")
    print(f"Calculated total genetic variance (Vg = H² * Vp): {H2} * {Vp} = {Vg}\n")

    print("--- Key Genetic Relationships ---")
    print("Total Genetic Variance (Vg) = Additive (Va) + Dominance (Vd) + Epistatic (Vi)")
    print("Narrow-sense heritability (h²) = Va / Vp")
    print("A standard Polygenic Score (PGS) primarily captures additive variance (Va), so the maximum variance it can explain is h².\n")

    # --- Analysis of Each Statement ---

    print("--- Evaluating Statement A ---")
    print("Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print("The maximum possible variance that ANY genetic predictor can explain is the total genetic variance, Vg.")
    print(f"The maximum proportion of variance explained is Vg / Vp = H² = {H2} (or {H2 * 100}%).")
    print("A PGS is a genetic predictor, so its explained variance must be less than or equal to H².")
    print("Therefore, a PGS cannot explain more than 50% of the variance.")
    print("CONCLUSION: Statement A is necessarily TRUE.\n")

    print("--- Evaluating Statement B ---")
    print("Statement B: Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print("This would only happen if h² approaches H², which means Va must be equal to Vg (i.e., Vd=0 and Vi=0).")
    print("Let's test a counterexample where non-additive effects exist:")
    # Scenario 1: Non-additive effects exist
    Va1 = 30.0
    Vd1 = 15.0
    Vi1 = 5.0
    print(f"Let Va = {Va1}, Vd = {Vd1}, Vi = {Vi1}. Note that {Va1} + {Vd1} + {Vi1} = {Vg}, which is consistent with H²=0.5.")
    h2_1 = Va1 / Vp
    print(f"In this case, the max variance a PGS can explain is h² = Va / Vp = {Va1} / {Vp} = {h2_1} (or {h2_1 * 100}%).")
    print(f"Since {h2_1 * 100}% is not 50%, the statement is not guaranteed.")
    print("CONCLUSION: Statement B is NOT necessarily true.\n")

    print("--- Evaluating Statement C ---")
    print("Statement C: Given an arbitrarily large GWAS, the polygenic score ... will not approach a variance explained of 50% due to ... non-linear effects.")
    print("This implies that non-additive effects (Vd or Vi) must exist. Let's test a counterexample where they don't.")
    # Scenario 2: All genetic variance is additive
    Va2 = 50.0
    Vd2 = 0.0
    Vi2 = 0.0
    print(f"Let Va = {Va2}, Vd = {Vd2}, Vi = {Vi2}. Note that {Va2} + {Vd2} + {Vi2} = {Vg}, which is consistent with H²=0.5.")
    h2_2 = Va2 / Vp
    print(f"In this case, h² = Va / Vp = {Va2} / {Vp} = {h2_2} (or {h2_2 * 100}%).")
    print("Here, it is theoretically possible for the PGS to approach explaining 50% of the variance.")
    print("CONCLUSION: Statement C is NOT necessarily true.\n")
    
    print("--- Evaluating Statement D ---")
    print("Statement D: The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("Epigenetic effects are typically considered part of the environmental variance (Ve).")
    print(f"The problem fixes the genetic contribution (Vg/Vp) at {H2}. Epigenetic effects contribute to the other 50%.")
    print("This does not place any new constraints on the components of Vg.")
    print(f"As shown in the analysis of Statement C, it is still possible for Va to equal Vg, making h² = H² = {H2}.")
    print("CONCLUSION: Statement D is NOT necessarily true.\n")

    print("--- Final Conclusion ---")
    print("Based on the analysis, only statement A is a necessary truth derived from the problem's premises.")

if __name__ == '__main__':
    analyze_heritability_statements()
    # Appending the final answer in the specified format for the platform.
    # The logic above concludes that A is the only necessarily true statement.
    sys.stdout.write("\n<<<A>>>\n")