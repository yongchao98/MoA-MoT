import sys

def reaction_analysis(eq_nBuLi=1.05):
    """
    Analyzes the stoichiometry of the borylation reaction and identifies side products.
    """
    # Let's assume a starting scale of 10.0 mmol for the aryl halide
    moles_aryl_halide = 0.010

    # Calculate moles of reagents based on stoichiometry
    moles_nBuLi_total = moles_aryl_halide * eq_nBuLi

    # The main reaction consumes n-BuLi in a 1:1 ratio with the aryl halide.
    moles_nBuLi_consumed = moles_aryl_halide
    moles_desired_product = moles_aryl_halide
    
    # Calculate the excess n-BuLi that causes the side reaction
    moles_nBuLi_excess = moles_nBuLi_total - moles_nBuLi_consumed
    
    # The excess n-BuLi forms the side product
    moles_side_product = moles_nBuLi_excess

    print("--- REACTION STOICHIOMETRY ANALYSIS ---")
    print(f"Using {eq_nBuLi:.2f} equivalents of n-BuLi:\n")
    print(f"1. Moles of starting material (Ar-I): {moles_aryl_halide:.4f} mol")
    print(f"2. Total moles of n-BuLi added: {moles_nBuLi_total:.4f} mol")
    print(f"3. Moles of n-BuLi consumed in main reaction: {moles_nBuLi_consumed:.4f} mol")
    print(f"4. Moles of n-BuLi remaining in excess: {moles_nBuLi_excess:.4f} mol\n")
    
    print("--- PRODUCTS FORMED ---")
    print("This leads to the formation of two different boron-containing products:")
    print(f" -> Desired Product ((2-bromo-4-chlorophenyl)boronic acid): {moles_desired_product:.4f} mol (Source of NMR Signal 1)")
    
    if moles_side_product > 1e-9: # Check if there is a meaningful amount of excess
        print(f" -> Side Product (n-butylboronic acid): {moles_side_product:.4f} mol (Source of NMR Signal 2)\n")
    else:
        print(" -> No significant side product is formed.\n")

    print("--- CHEMICAL EQUATIONS ---")
    print("The numbers in the equations represent the stoichiometric coefficients.")
    print("\nMain Reaction Pathway:")
    print("1 Ar-I + 1 n-BuLi --> 1 Ar-Li + 1 n-BuI")
    print("1 Ar-Li + 1 B(OMe)3 --(Workup)--> 1 Ar-B(OH)2\n")

    print("Side Reaction Pathway (due to excess n-BuLi):")
    print("1 n-BuLi + 1 B(OMe)3 --(Workup)--> 1 n-Bu-B(OH)2\n")

    print("--- CONCLUSION ---")
    print("The presence of excess n-BuLi leads to a side product, explaining the two B NMR signals.")
    print("To solve this, a more precise amount of n-BuLi should be used (Option C).")
    
# Run the analysis for the described problem
reaction_analysis(eq_nBuLi=1.05)

<<<C>>>