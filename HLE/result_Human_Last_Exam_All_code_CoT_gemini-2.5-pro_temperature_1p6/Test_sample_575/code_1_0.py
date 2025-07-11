import sys

def heritability_analysis():
    """
    Analyzes statements about heritability and polygenic scores to find the one
    that is necessarily true given the problem's premises.
    """

    # --- Problem Setup & Definitions ---
    H_sq = 0.5  # Broad-sense heritability (given)

    print("--- Analysis of Heritability Statements ---")
    print(f"The problem states that the broad-sense heritability (H^2) of a phenotype is {H_sq}.")
    print("\nKey Definitions:")
    print("1. Broad-sense Heritability (H^2) = Vg / Vp = (Va + Vd + Vi) / Vp")
    print("   (Vg=Total Genetic Variance, Va=Additive, Vd=Dominance, Vi=Epistasis, Vp=Phenotypic Variance)")
    print("2. Narrow-sense Heritability (h^2) = Va / Vp")
    print("3. Polygenic Score (PGS) from a standard GWAS primarily captures additive effects.")
    print("   Therefore, the variance a PGS can explain (R^2) approaches h^2 in an infinitely large study.")
    print(f"\nFundamental Constraint 1: By definition, h^2 <= H^2. So, h^2 must be <= {H_sq}.")
    print(f"Fundamental Constraint 2: The absolute maximum variance any genetic predictor can explain is H^2, which is {H_sq*100}%.")
    print("-" * 40)

    # --- Modeling Scenarios ---
    print("\nLet's test the statements against two plausible scenarios:\n")

    # Scenario 1: A purely additive trait.
    print("--- Scenario 1: Genetic variance is purely additive (Vd=0, Vi=0) ---")
    # In this case, Vg = Va, which means h^2 = H^2.
    h_sq_scenario1 = H_sq
    pgs_r_sq_1 = h_sq_scenario1
    print(f"   In this case, the narrow-sense heritability (h^2) is equal to H^2.")
    print(f"   h^2 = {h_sq_scenario1}")
    print(f"   The variance explained by an ideal PGS will approach {pgs_r_sq_1*100}%.\n")


    # Scenario 2: A trait with non-additive effects.
    print("--- Scenario 2: Genetic variance has non-additive components (Vd>0 or Vi>0) ---")
    # Let's assume additive variance is 60% of total genetic variance.
    additive_fraction = 0.6
    # In this case, Va < Vg, which means h^2 < H^2.
    h_sq_scenario2 = H_sq * additive_fraction
    pgs_r_sq_2 = h_sq_scenario2
    print(f"   Let's assume Va is {additive_fraction*100}% of Vg. The remaining {round((1-additive_fraction)*100)}% is from dominance/epistasis.")
    print(f"   In this case, the narrow-sense heritability (h^2) = {H_sq} * {additive_fraction} = {h_sq_scenario2}")
    print(f"   The variance explained by an ideal PGS will approach {pgs_r_sq_2*100}%.\n")
    print("-" * 40)


    # --- Evaluating the Answer Choices ---
    print("\nEvaluating each statement based on the scenarios:\n")

    # Statement A
    print("A. The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print(f"   - In Scenario 1, PGS explains {pgs_r_sq_1*100}%. This is not more than 50%.")
    print(f"   - In Scenario 2, PGS explains {pgs_r_sq_2*100}%. This is not more than 50%.")
    print("   - This holds true because the total heritable variance (H^2) is 50%, which is the absolute ceiling for any genetic predictor.")
    print("   - Verdict: NECESSARILY TRUE.\n")

    # Statement B
    print("B. Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print(f"   - This is TRUE for Scenario 1 (PGS approaches {pgs_r_sq_1*100}%).")
    print(f"   - This is FALSE for Scenario 2 (PGS approaches {pgs_r_sq_2*100}%).")
    print("   - Since the statement can be false, it is not necessarily true.")
    print("   - Verdict: NOT NECESSARILY TRUE.\n")

    # Statement C
    print("C. ...the polygenic score...will not approach...50% due to...non-linear effects...")
    print(f"   - This is FALSE for Scenario 1, where there are no non-linear effects and the PGS does approach 50%.")
    print(f"   - This is TRUE for Scenario 2.")
    print("   - Since the statement can be false, it is not necessarily true.")
    print("   - Verdict: NOT NECESSARILY TRUE.\n")

    # Statement D
    print("D. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("   - h^2 is always <= H^2 (0.5) by definition. The question is if epigenetics forces it to be strictly less than 0.5.")
    print("   - This is false. A trait could be purely additive (h^2 = 0.5) while still being influenced by epigenetics (which would contribute to non-genetic variance).")
    print("   - Verdict: NOT NECESSARILY TRUE.\n")

# Run the analysis
heritability_analysis()
<<<A>>>