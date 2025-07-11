import sys
# Redirect print to a string buffer to control final output format.
# The user sees the final printed output, not this setup code.
import io
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def analyze_heritability():
    """
    This function models the components of heritability to evaluate the statements.
    """
    # --- 1. Given Information & Definitions ---
    # We are given the broad-sense heritability (H^2).
    H2 = 0.5

    # For demonstration, let's set a total phenotypic variance (Vp).
    # The actual value doesn't matter, as heritability is a ratio.
    Vp = 100

    # Total Genetic Variance (Vg) is defined by H^2.
    Vg = H2 * Vp

    # Environmental Variance (Ve) is the remaining variance.
    Ve = Vp - Vg

    print("--- Scenario Setup Based on H^2 = 0.5 ---")
    print(f"Broad-sense Heritability (H^2): {H2}")
    print(f"Total Genetic Variance (Vg): {Vg} (which is 50.0% of Vp)")
    print(f"Total Environmental Variance (Ve): {Ve} (which is 50.0% of Vp)")
    print("-" * 20)

    # --- 2. A Plausible Partition of Genetic Variance for a Polygenic Trait ---
    # For a complex trait, genetic variance (Vg) is composed of:
    # Va (Additive), Vd (Dominance), and Vi (Epistatic/Interaction).
    # Typically, Va < Vg. Let's model this.
    # We'll assume 70% of genetic variance is additive.
    Va = 0.7 * Vg
    Vd = 0.2 * Vg
    Vi = 0.1 * Vg

    # Narrow-sense heritability (h^2) is the proportion of variance from additive effects.
    h2 = Va / Vp

    # A standard Polygenic Score (PGS) from GWAS is a linear model that primarily
    # captures additive effects. Its predictive power (R^2) is limited by h^2.
    max_pgs_r_squared = h2

    print("--- Modeling Components and PGS Performance ---")
    print(f"Total Genetic Variance (Vg) is partitioned into:")
    print(f"  - Additive Variance (Va): {Va}")
    print(f"  - Dominance Variance (Vd): {Vd}")
    print(f"  - Epistatic Variance (Vi): {Vi}")
    print(f"Check sum: Va + Vd + Vi = {Va + Vd + Vi}")
    print("\nThis leads to:")
    print(f"Narrow-sense Heritability (h^2 = Va / Vp): {h2}")
    print(f"Max R^2 of a linear PGS is limited by h^2, so it approaches: {max_pgs_r_squared:.2f}")
    print("-" * 20)

    # --- 3. Evaluating Each Statement ---
    print("--- Evaluation of Answer Choices ---")

    # Statement A: PGS R^2 cannot be > 50%
    is_A_true = max_pgs_r_squared <= H2
    print(f"A. The PGS cannot explain more than {H2*100}% of variance.")
    print(f"   Logic: The PGS is limited by h^2 ({h2:.2f}), and h^2 is always <= H^2 ({H2:.2f}).")
    print(f"   Conclusion: Statement A is necessarily TRUE.\n")

    # Statement B: PGS R^2 will approach 50%
    # This implies max_pgs_r_squared would approach H2
    is_B_true = h2 == H2
    print(f"B. The PGS will approach a variance explained of {H2*100}%.")
    print(f"   Logic: The PGS approaches h^2 ({h2:.2f}), not H^2 ({H2:.2f}). This is only true if non-additive variance is zero.")
    print(f"   Conclusion: Statement B is FALSE.\n")

    # Statement C: PGS R^2 will not approach 50% due to non-linear effects
    is_C_true = (Vd + Vi) > 0
    print(f"C. The PGS will not approach {H2*100}% due to dominance/epistasis.")
    print(f"   Logic: Because non-additive variance (Vd+Vi = {Vd+Vi}) is > 0, h^2 ({h2:.2f}) is strictly less than H^2 ({H2:.2f}).")
    print(f"   Conclusion: Statement C is TRUE.\n")

    # Statement D: Epigenetics would limit h^2 to < 50%
    print(f"D. Epigenetic effects would limit h^2 to be < {H2}.")
    print(f"   Logic: Epigenetic effects contribute to Environmental Variance (Ve), not Genetic Variance (Vg).")
    print(f"   The value of Ve does not determine the internal partitioning of Vg into Va, Vd, Vi.")
    print(f"   Conclusion: Statement D is FALSE.\n")

    print("--- Final Result ---")
    print("Statements A and C are the correct conclusions from the premises.")

# Execute the analysis
analyze_heritability()
# Restore standard output and print the captured analysis
sys.stdout = old_stdout
print(captured_output.getvalue())