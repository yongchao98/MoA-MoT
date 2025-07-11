import sys

def analyze_heritability_statements():
    """
    Analyzes statements about heritability and polygenic scores based on given principles.
    This script serves to lay out the logical argument step-by-step.
    """

    # --- Setup from the problem statement ---
    # Broad-sense heritability (H2) is the proportion of phenotypic variance (Vp)
    # explained by total genetic variance (Vg).
    # H2 = Vg / Vp
    H2 = 0.5

    print("--- Principles of Quantitative Genetics ---")
    print("Vp: Total Phenotypic Variance")
    print("Vg: Total Genetic Variance")
    print("Va: Additive Genetic Variance (effects that sum up)")
    print("Vd: Dominance Variance (non-linear interactions at one locus)")
    print("Vi: Epistatic Variance (non-linear interactions between loci)")
    print("\nKey Equations:")
    print("1. Vg = Va + Vd + Vi")
    print(f"2. Broad-Sense Heritability (H2) = Vg / Vp = {H2} (given)")
    print("3. Narrow-Sense Heritability (h2) = Va / Vp")
    print("4. A standard linear Polygenic Score (PGS) explains variance approximately equal to h2.")
    print("-" * 40)

    # --- Deriving the core inequality ---
    print("--- Deriving a Core Relationship ---")
    print("From 'Vg = Va + Vd + Vi', we know that Va <= Vg (since Vd and Vi cannot be negative).")
    print("Dividing by Vp: (Va / Vp) <= (Vg / Vp)")
    print("This means: h2 <= H2")
    print(f"Given H2 = {H2}, it must be true that h2 <= {H2}.")
    print("Since the maximum variance a linear PGS can explain is h2, let's call it V_PGS.")
    print(f"Therefore, a fundamental truth is: V_PGS <= h2 <= {H2}.")
    print("-" * 40)

    # --- Evaluating the Statements ---
    print("--- Evaluating Answer Choices ---")

    # Statement A
    print("\n[A] The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print("This statement translates to: V_PGS <= 0.5?")
    print(f"From our derivation, we established that V_PGS <= h2 <= H2, and H2 = {H2}.")
    print(f"So, it is necessarily true that V_PGS <= {H2}. The variance explained by the genetic score cannot exceed the total variance explained by all genetic factors.")
    is_A_true = True
    print("Conclusion: Statement A is TRUE.")

    # Statement B
    print("\n[B] Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%.")
    print("This statement translates to: V_PGS will approach 0.5?")
    print("This would only be true if h2 = H2 = 0.5.")
    print("For h2 to equal H2, additive variance (Va) must equal total genetic variance (Vg).")
    print("This would require non-additive variance (Vd + Vi) to be exactly 0.")
    print("We cannot assume that non-additive variance is zero for a polygenic trait. It's theoretically possible, but not necessarily true.")
    is_B_true = False
    print("Conclusion: Statement B is NOT necessarily true.")

    # Statement C
    print("\n[C] Given an arbitrarily large GWAS, the polygenic score constructed via linearly summing GWAS effect sizes will not approach a variance explained of 50% due to gene-gene interactions and other non-linear effects such as dominance.")
    print("This statement claims V_PGS will approach a value < 0.5, because the linear model of a PGS cannot capture non-linear effects (Vd and Vi).")
    print("This is the flip side of statement B. It asserts that non-additive variance (Vd + Vi) is greater than 0.")
    print("For any realistic polygenic trait, it is biologically expected that there are some non-additive effects, meaning Vd + Vi > 0.")
    print("If Vd + Vi > 0, then Va < Vg, which means h2 < H2.")
    print(f"In this case, V_PGS (which approaches h2) would be less than H2 ({H2}).")
    print("This statement correctly identifies the primary limitation of standard PGS in capturing the full genetic contribution.")
    is_C_true = True
    print("Conclusion: Statement C is considered TRUE in the context of real-world polygenic architecture.")

    # Statement D
    print("\n[D] The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5.")
    print("We already know that h2 <= 0.5 by definition. The statement claims a strict inequality, h2 < 0.5, is forced by epigenetic effects.")
    print("Heritable epigenetic effects could be complex, but if we model non-heritable epigenetic variance as part of the Environmental Variance (Ve), it increases the denominator (Vp) of h2 = Va/Vp, thus lowering h2.")
    print("However, it is theoretically possible to construct a scenario where epigenetic effects exist, but h2 = H2 = 0.5 (if Vd=Vi=0). The problem gives us the final state where H2=0.5, not what happens when effects are added.")
    print("Therefore, the presence of epigenetic effects does not *necessarily* force h2 to be strictly less than 0.5.")
    is_D_true = False
    print("Conclusion: Statement D is NOT necessarily true.")
    print("-" * 40)

    # --- Final Summary ---
    print("\n--- Summary ---")
    print(f"Statement A is necessarily true by definition.")
    print(f"Statement C is true because it correctly describes a fundamental limitation of the specified method (linear PGS) when applied to complex polygenic traits, which are assumed to have non-additive effects.")
    print("Statements B and D are not necessarily true.")
    print("Therefore, the correct choice is that only A and C are necessarily true.")

if __name__ == '__main__':
    analyze_heritability_statements()