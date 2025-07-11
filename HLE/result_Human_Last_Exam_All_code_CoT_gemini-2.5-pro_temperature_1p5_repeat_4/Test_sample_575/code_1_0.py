import sys
import io

# Redirect stdout to capture the output for the final answer block
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def analyze_heritability():
    """
    Analyzes the relationship between broad-sense heritability (H^2),
    narrow-sense heritability (h^2), and the variance explained by a polygenic score (PGS).
    """
    # --- 1. Define Parameters from the Problem ---
    H2 = 0.5  # Given broad-sense heritability

    # For demonstration, we assume a total phenotypic variance (Vp).
    # The specific value doesn't matter as we are dealing with proportions,
    # but using 100 makes percentages intuitive.
    Vp = 100.0

    # --- 2. Calculate Total Genetic Variance (Vg) ---
    # From the definition H^2 = Vg / Vp
    Vg = H2 * Vp

    print("--- Core Concepts & Given Information ---")
    print(f"Broad-Sense Heritability (H^2): {H2}")
    print(f"This means that all genetic factors (Vg) combined explain {H2*100}% of the total phenotypic variance (Vp).")
    print(f"Assuming Vp = {Vp}, the Total Genetic Variance (Vg) = {H2} * {Vp} = {Vg}")
    print("The total genetic variance (Vg) is composed of Vg = Va + Vd + Vi")
    print("(Va=Additive, Vd=Dominance, Vi=Epistasis)")
    print("\nA Polygenic Score (PGS) from GWAS primarily captures additive effects (Va).")
    print("The maximum variance a perfect PGS can explain is the narrow-sense heritability (h^2 = Va / Vp).")
    print("-" * 50)

    # --- 3. Model Two Scenarios to Test the Statements ---

    # Scenario 1: A purely additive trait
    # Here, all genetic variance is additive. Vd and Vi are zero.
    Va_1 = Vg
    Vd_1 = 0.0
    Vi_1 = 0.0
    h2_1 = Va_1 / Vp

    print("Scenario 1: Genetic architecture is PURELY ADDITIVE")
    print(f"  Vg = Va + Vd + Vi")
    print(f"  {Vg} = {Va_1} + {Vd_1} + {Vi_1}")
    print(f"  Narrow-sense heritability (h^2) = Va / Vp = {Va_1} / {Vp} = {h2_1}")
    print(f"  Max variance explained by PGS in this case: {h2_1*100}%")
    print("-" * 50)


    # Scenario 2: A trait with non-additive effects
    # Here, Vg is split between additive, dominance, and epistasis.
    # The split is arbitrary, as long as Va + Vd + Vi = Vg.
    Va_2 = 30.0
    Vd_2 = 15.0
    Vi_2 = 5.0 # Check: 30 + 15 + 5 = 50 = Vg
    h2_2 = Va_2 / Vp

    print("Scenario 2: Genetic architecture includes NON-ADDITIVE effects")
    print(f"  Vg = Va + Vd + Vi")
    print(f"  {Vg} = {Va_2} + {Vd_2} + {Vi_2}")
    print(f"  Narrow-sense heritability (h^2) = Va / Vp = {Va_2} / {Vp} = {h2_2}")
    print(f"  Max variance explained by PGS in this case: {h2_2*100}%")
    print("-" * 50)

    # --- 4. Evaluate Each Statement Based on the Scenarios ---
    print("\n--- Evaluating Answer Choices ---")
    print("The key relationship is: PGS_explained_variance <= h^2 <= H^2")
    print(f"Given H^2 = 0.5, we know PGS_explained_variance <= {H2}.\n")

    # A. The polygenic score can not explain more than 50% of the variance in the phenotype.
    print("A: This states PGS_explained_variance <= 0.5. This is ALWAYS TRUE because h^2 <= H^2 (0.5).")
    print(f"   In Scenario 1, PGS explains {h2_1*100}%. In Scenario 2, PGS explains {h2_2*100}%. Neither is > 50%.")
    print("   This statement is NECESSARILY TRUE.\n")

    # B. Given an arbitrarily large GWAS, the polygenic score will approach a variance explained of 50%
    print("B: This states the explained variance will approach 0.5. This requires h^2 = 0.5.")
    print("   This is true in Scenario 1, but false in Scenario 2. So, it is NOT NECESSARILY TRUE.\n")

    # C. ...the polygenic score... will not approach a variance explained of 50%...
    print("C: This states the explained variance will be < 0.5. This requires h^2 < 0.5.")
    print("   This is true in Scenario 2, but false in Scenario 1. So, it is NOT NECESSARILY TRUE.\n")

    # D. The existence of any epigenetic effects would limit the narrow-sense heritability to less than 0.5
    print("D: This is incorrect. The problem defines H^2 = 0.5 regardless of its source. The logical relationship")
    print("   h^2 <= H^2 still holds. It does not force h^2 to be strictly less than H^2. So, NOT NECESSARILY TRUE.\n")

    # E. None of the other answer choices are correct
    print("E: This is incorrect because statement A is necessarily true.")


# Run the analysis
analyze_heritability()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())