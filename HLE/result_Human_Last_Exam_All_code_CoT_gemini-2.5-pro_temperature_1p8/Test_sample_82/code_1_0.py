import pandas as pd

def analyze_heritability():
    """
    Analyzes the statements about heritability and polygenic scores
    based on the provided information.
    """
    # --- Initial Setup based on the problem ---
    H2 = 0.5  # Broad-sense heritability
    Vp = 1.0  # Assume total phenotypic variance is 1 for simplicity
    Vg = H2 * Vp  # Total genetic variance

    print("--- Problem Setup ---")
    print(f"Broad-sense heritability (H2): {H2}")
    print(f"Total phenotypic variance (Vp): {Vp:.1f}")
    print(f"Total genetic variance (Vg = H2 * Vp): {Vg:.1f}")
    print(f"Environmental variance (Ve = Vp - Vg): {1.0 - Vg:.1f}")
    print("-" * 25)

    # --- Scenario 1: All genetic variance is additive (special case) ---
    print("--- Scenario 1: Purely Additive Genetic Variance ---")
    s1_Va = Vg
    s1_Vd_Vi = 0.0
    s1_h2 = s1_Va / Vp
    print(f"In this scenario, we assume Vg = Va.")
    print(f"Additive variance (Va): {s1_Va:.1f}")
    print(f"Non-additive variance (Vd + Vi): {s1_Vd_Vi:.1f}")
    print(f"Narrow-sense heritability (h2 = Va / Vp): {s1_h2:.1f}")
    print(f"Max R-squared a linear PGS can explain: {s1_h2:.1f} (approaches 50%)")
    print("-" * 25)

    # --- Scenario 2: Genetic variance has non-additive components (typical case for polygenic traits) ---
    print("--- Scenario 2: Mixed Additive & Non-Additive Variance ---")
    # We assign some plausible values where Va < Vg
    s2_Va = 0.3
    s2_Vd_Vi = Vg - s2_Va
    s2_h2 = s2_Va / Vp
    print(f"In this scenario, common for polygenic traits, Vg = Va + Vd + Vi.")
    print(f"Additive variance (Va): {s2_Va:.1f}")
    print(f"Non-additive variance (Vd + Vi): {s2_Vd_Vi:.2f}")
    print(f"Narrow-sense heritability (h2 = Va / Vp): {s2_h2:.1f}")
    print(f"Max R-squared a linear PGS can explain: {s2_h2:.1f} (approaches 30%)")
    print("-" * 25)

    # --- Evaluating the statements ---
    print("--- Final Evaluation of Statements ---")
    # Statement A
    print("A. The polygenic score can not explain more than 50% of the variance.")
    print("   - Evaluation: The total genetic variance (Vg) is 50% (H2=0.5). A genetic predictor cannot explain")
    print("     more variance than exists. The R-squared of a PGS is always <= H2.")
    print("   - Verdict: Necessarily TRUE.\n")

    # Statement B
    print("B. Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print("   - Evaluation: This is only true in Scenario 1, where all genetic variance is additive (h2 = H2).")
    print("     It is false in Scenario 2. Since we cannot assume Scenario 1, this is not *necessarily* true.")
    print("   - Verdict: Not necessarily true.\n")
    
    # Statement C
    print("C. Given an arbitrarily large GWAS, the polygenic score ... will not approach ... 50% due to ... non-linear effects.")
    print("   - Evaluation: A linear PGS is limited by narrow-sense heritability (h2). The term 'polygenic' implies")
    print("     a complex trait where non-additive effects (dominance, epistasis) exist, making h2 < H2.")
    print("     Therefore, the PGS will be limited to a value less than 50%.")
    print("   - Verdict: Necessarily TRUE (in the context of a polygenic trait).\n")
    
    # Statement D
    print("D. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("   - Evaluation: Epigenetic effects contribute to environmental variance (Ve), not genetic variance (Vg).")
    print("     The composition of Vg is independent of Ve. h2 could still equal H2 (0.5).")
    print("   - Verdict: Not necessarily true.\n")
    
    print("--- Conclusion ---")
    print("Statements A and C are the only ones that are necessarily true.")

analyze_heritability()
<<<E>>>