def analyze_boronic_acid_synthesis(eq_borate):
    """
    Calculates the amount of product and byproduct for a boronic acid synthesis.

    Args:
        eq_borate (float): The equivalents of trimethyl borate used relative to the starting material.
    """
    # Assume we start with 1.0 mole of the limiting reagent for simplicity
    moles_starting_material = 1.0
    
    # The reaction stoichiometry for product formation is 1:1
    moles_product_formed = moles_starting_material
    moles_borate_consumed = moles_product_formed
    
    # Calculate total and excess borate used
    moles_borate_total = moles_starting_material * eq_borate
    moles_borate_excess = moles_borate_total - moles_borate_consumed
    
    # Excess trimethyl borate is hydrolyzed to boric acid byproduct during workup
    # B(OMe)3 (excess) + 3 H2O -> B(OH)3 + 3 MeOH
    moles_boric_acid_byproduct = moles_borate_excess
    
    print(f"Analysis for {eq_borate} eq of Trimethyl Borate:")
    print(f"  - Starting Material: {moles_starting_material:.1f} mol")
    print(f"  - Desired Product (ArB(OH)2) formed: {moles_product_formed:.1f} mol")
    
    # Print the final equation for byproduct formation with calculated numbers
    print("  - Byproduct formation equation upon workup:")
    print(f"    {moles_borate_excess:.1f} mol B(OMe)3 (excess) + H2O -> {moles_boric_acid_byproduct:.1f} mol B(OH)3 (byproduct)")
    
    print(f"  - Molar ratio of Byproduct (B(OH)3) to Product (ArB(OH)2): {moles_boric_acid_byproduct:.1f} : {moles_product_formed:.1f}\n")

# --- Main ---
print("The problem of two B-NMR signals is likely due to a large amount of boric acid byproduct.")
print("This byproduct comes from the hydrolysis of the large excess of trimethyl borate used.\n")

# Scenario 1: Original Procedure
print("--- SCENARIO 1: ORIGINAL PROCEDURE ---")
analyze_boronic_acid_synthesis(eq_borate=5.0)

# Scenario 2: Proposed Solution
print("--- SCENARIO 2: PROPOSED SOLUTION (Option D) ---")
analyze_boronic_acid_synthesis(eq_borate=1.5)

print("CONCLUSION:")
print("Using a large excess (5 eq) of trimethyl borate results in 4 moles of boric acid byproduct for every 1 mole of product.")
print("This makes purification difficult and leads to multiple NMR signals.")
print("Using less trimethyl borate (e.g., 1.5 eq) drastically reduces the byproduct to only 0.5 moles, leading to a much cleaner product.")
