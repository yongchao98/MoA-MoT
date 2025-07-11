import sys

def analyze_genetics_problem():
    """
    Analyzes the heritability problem statement by statement.
    """
    # --- Step 1: Define variables from the problem ---
    # H2 is the broad-sense heritability. It represents the total proportion of
    # phenotypic variance (Vp) that is due to all genetic variance (Vg).
    # Vg includes additive (Va), dominance (Vd), and epistatic (Vi) effects.
    # H2 = Vg / Vp
    H2 = 0.5

    print("--- Problem Setup ---")
    print(f"Broad-sense heritability (H²) is given as: {H2}")
    print("This means the total genetic variance (Vg) accounts for 50% of the total phenotypic variance (Vp).")
    print("The remaining 50% is environmental variance (Ve).\n")
    print("The key relationships are:")
    print("  - Vp = Vg + Ve")
    print("  - Vg = Va + Vd + Vi (Additive + Dominance + Epistatic variance)")
    print(f"  - H² = Vg / Vp = {H2}")
    print("  - h² = Va / Vp (Narrow-sense heritability)")
    print("  - A standard linear Polygenic Score (PGS) predicts phenotype using additive effects, so its maximum explanatory power approaches h².\n")


    # --- Step 2: Evaluate each statement ---
    print("--- Evaluating the Statements ---")

    # Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype.
    # The PGS is a predictor based on genetic information.
    # The total variance attributable to ALL genetic factors (Vg) is H² = 0.5 or 50%.
    # Therefore, no model based purely on genetics can explain more than 50% of the variance.
    # This statement must be true.
    print("Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype.")
    is_A_true = True
    print(f"Reasoning: The maximum variance explainable by any genetic predictor is the total genetic variance, which is H² = {H2}. A PGS is a genetic predictor, so it cannot explain more than {H2*100}% of the variance.")
    print(f"Conclusion: Statement A is necessarily TRUE.\n")


    # Statement B: Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.
    # A standard PGS sums linear, additive effects. Its predictive power approaches h² = Va/Vp.
    # The total genetic variance is H² = Vg/Vp = 0.5.
    # The PGS would only approach 0.5 if h² = H², which means Vd and Vi are both 0.
    # The problem does not state that non-additive effects are absent.
    # Therefore, we cannot be certain it will approach 0.5. It will approach h², which could be less than 0.5.
    print("Statement B: Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    is_B_true = False
    print(f"Reasoning: A standard PGS captures additive variance (Va), so its accuracy approaches h² (Va/Vp). It only approaches H² ({H2}) if all genetic variance is additive (h²=H²). This is not guaranteed.")
    print(f"Conclusion: Statement B is NOT necessarily true.\n")


    # Statement C: Given an arbitrarily large GWAS, the polygenic score constructed via linearly summing
    # GWAS effect sizes will not approach a variance explained of 50% due to gene-gene interactions
    # and other non-linear effects such as dominance.
    # This is the converse of B. It asserts that because of non-additive effects (Vd, Vi), Va < Vg.
    # If Va < Vg, then h² < H². The PGS accuracy approaches h², so it would not reach H² (0.5).
    # For a complex polygenic trait, the existence of non-additive effects is a standard assumption.
    # Without them, genetics would be purely additive, which is biologically rare for complex traits.
    print("Statement C: Given an arbitrarily large GWAS, the polygenic score constructed via linearly summing GWAS effect sizes will not approach a variance explained of 50% due to gene-gene interactions and other non-linear effects such as dominance.")
    is_C_true = True
    print(f"Reasoning: This statement assumes that non-additive genetic effects (dominance, epistasis) exist for a polygenic trait, which is a standard and biologically realistic assumption. If they exist, then additive variance (Va) is less than total genetic variance (Vg). This means h² < H². The PGS accuracy approaches h², so it will not reach the H² value of {H2}.")
    print(f"Conclusion: Statement C is considered necessarily TRUE under standard biological assumptions.\n")


    # Statement D: The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.
    # Epigenetic effects that are not transgenerationally inherited are typically considered part of the environmental variance (Ve).
    # The value H²=0.5 is given. This means Vg/Vp=0.5.
    # This implies h² = Va/Vp <= Vg/Vp = 0.5. So h² is already limited to a maximum of 0.5.
    # The existence of epigenetic effects (in Ve) does not force h² to be strictly *less than* 0.5.
    # It's possible for all genetic variance to be additive (h²=0.5) while epigenetic effects still exist (as part of Ve).
    print("Statement D: The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    is_D_true = False
    print(f"Reasoning: Narrow-sense heritability (h²) is already limited to be less than or equal to H² ({H2}). Epigenetic effects are generally considered part of environmental variance (Ve). The existence of Ve doesn't force h² to be strictly *less than* {H2}. A scenario where h² = {H2} and Ve > 0 (containing epigenetic effects) is possible.")
    print(f"Conclusion: Statement D is NOT necessarily true.\n")

    # --- Step 3: Final Conclusion ---
    print("--- Final Conclusion ---")
    print("Based on the analysis:")
    print(f"- Statement A is TRUE.")
    print(f"- Statement B is FALSE.")
    print(f"- Statement C is TRUE.")
    print(f"- Statement D is FALSE.")
    print("\nTherefore, the correct choice is the one that states only A and C are correct.")

if __name__ == "__main__":
    analyze_genetics_problem()
<<<E>>>