import sys
import io

# A helper function to redirect print to a string
def execute_and_capture_output(func):
    old_stdout = sys.stdout
    sys.stdout = captured_output = io.StringIO()
    func()
    sys.stdout = old_stdout
    return captured_output.getvalue()

def solve():
    """
    Analyzes the quantitative genetics problem and demonstrates the logic.
    """
    # --- Setup from the problem ---
    # In a theoretically ideal population, a phenotype has a broad-sense heritability of 0.5.
    H2 = 0.5  # Broad-sense heritability

    # For the purpose of demonstration, let's assume a total phenotypic variance (Vp).
    # The actual value doesn't matter, so we'll use 100 for easy interpretation of percentages.
    Vp = 100

    print("Step 1: Define the given parameters based on the problem statement.")
    print(f"Given Broad-sense heritability (H2) = {H2}")
    print(f"Let's assume a total phenotypic variance (Vp) for our model = {Vp}\n")

    # From H2 = Vg / Vp, we can calculate the total genetic variance (Vg)
    Vg = H2 * Vp
    print("Step 2: Calculate the total genetic variance (Vg).")
    print("Formula: Vg = H2 * Vp")
    print(f"Calculation: Vg = {H2} * {Vp} = {Vg}\n")

    # --- Analysis of Statement A ---
    print("--- Analysis of Statement A ---")
    print("Statement A: The polygenic score can not explain more than 50% of the variance in the phenotype.")
    print("A polygenic score (PGS) is built from genetic data.")
    print("The maximum possible phenotypic variance that can be explained by any genetic factor is the total genetic variance, Vg.")
    print(f"The maximum proportion of variance explained by genetics is Vg / Vp, which is the definition of H2.")
    print(f"Therefore, the variance explained by the PGS must be <= H2, which is {H2}.")
    print("Conclusion: Statement A is TRUE.\n")

    # --- Analysis of Statement C ---
    print("--- Analysis of Statement C ---")
    print("Statement C: Given an arbitrarily large GWAS, the polygenic score constructed via linearly summing GWAS effect sizes will not approach a variance explained of 50% due to gene-gene interactions and other non-linear effects such as dominance.")
    print("This statement addresses the difference between broad-sense (H2) and narrow-sense (h2) heritability.")
    print("Total genetic variance (Vg) is composed of additive (Va) and non-additive (Vd, Vi) parts: Vg = Va + Vd + Vi.")
    print("A standard linear PGS estimates h2, which is calculated as: h2 = Va / Vp.")
    print("Let's model a general scenario where the genetic variance is not purely additive.")
    
    # We assign plausible values to the components of Vg.
    # We know Vg must sum to 50 in our model.
    # Let's say Va = 30 and the non-additive part (Vd + Vi) = 20.
    Va_scenario = 30
    Vd_Vi_scenario = 20
    Vg_scenario_check = Va_scenario + Vd_Vi_scenario
    
    print(f"\nLet's assume additive variance (Va) = {Va_scenario} and non-additive variance (Vd + Vi) = {Vd_Vi_scenario}.")
    print(f"(Check: Va + Vd + Vi = {Va_scenario} + {Vd_Vi_scenario} = {Vg_scenario_check}, which matches our calculated Vg of {Vg}).")
    
    # Now, calculate the narrow-sense heritability (h2) for this scenario.
    h2_scenario = Va_scenario / Vp
    
    print("\nIn this scenario, the variance explained by the PGS (which approaches h2) would be:")
    print("Formula: h2 = Va / Vp")
    print(f"Calculation: h2 = {Va_scenario} / {Vp} = {h2_scenario}")
    print(f"\nThe PGS would approach explaining {h2_scenario * 100}% of the variance, not the full {H2 * 100}% (H2).")
    print("The reason for this shortfall is the non-additive variance, which the linear PGS model does not capture.")
    print("Since this limitation is inherent to the method for any trait with non-additive effects (the general case), Statement C is also a necessary truth about the methodology.")
    print("Conclusion: Statement C is TRUE.\n")
    
    print("--- Final Conclusion ---")
    print("Both statements A and C are necessarily true.")
    print("Statement B is incorrect because h2 is generally less than H2.")
    print("Statement D is incorrect because epigenetic effects do not force h2 to be less than H2.")
    print("The correct choice is the one that includes both A and C.")

# Execute the demonstration and print the output.
output = execute_and_capture_output(solve)
print(output)
<<<E>>>