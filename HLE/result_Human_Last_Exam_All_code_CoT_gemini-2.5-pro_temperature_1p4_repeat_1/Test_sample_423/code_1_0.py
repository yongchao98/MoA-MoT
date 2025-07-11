def solve_structure_x():
    """
    Solves the chemical puzzle to determine the structure of substance X.
    The function prints the step-by-step reasoning and calculations.
    """

    # --- Step 1: Combustion Analysis of X ---
    print("--- Step 1: Combustion Analysis of X ---")
    mass_co2 = 0.7472  # g
    mass_h2o = 0.1834  # g
    molar_mass_co2 = 12.011 + 2 * 15.999
    molar_mass_h2o = 2 * 1.008 + 15.999

    moles_co2 = mass_co2 / molar_mass_co2
    moles_c_in_x = moles_co2
    mass_c_in_x = moles_c_in_x * 12.011

    moles_h2o = mass_h2o / molar_mass_h2o
    moles_h_in_x = 2 * moles_h2o
    mass_h_in_x = moles_h_in_x * 1.008

    print(f"Moles of C in sample = {moles_c_in_x:.5f} mol")
    print(f"Moles of H in sample = {moles_h_in_x:.5f} mol")

    # Ratio of H to C
    h_c_ratio = moles_h_in_x / moles_c_in_x
    print(f"Ratio of H:C = {h_c_ratio:.3f}:1.0")
    print("This ratio is very close to 1.2:1, which is equivalent to 6:5 or 12:10.")
    print("This suggests an empirical formula component of (C5H6)n or (C10H12)n.")
    print("\n")


    # --- Step 2: Analysis of Substance B ---
    print("--- Step 2: Analysis of Substance B ---")
    mass_frac_c_b = 0.50
    mass_frac_h_b = 0.10
    mass_frac_o_b = 0.40

    # Relative moles in 100g of B
    moles_c_b = mass_frac_c_b / 12.011
    moles_h_b = mass_frac_h_b / 1.008
    moles_o_b = mass_frac_o_b / 15.999

    # Find the simplest ratio
    min_moles_b = min(moles_c_b, moles_h_b, moles_o_b)
    c_ratio_b = moles_c_b / min_moles_b
    h_ratio_b = moles_h_b / min_moles_b
    o_ratio_b = moles_o_b / min_moles_b

    print(f"Raw ratio C:H:O = {c_ratio_b:.2f}:{h_ratio_b:.2f}:{o_ratio_b:.2f}")
    # Multiply by 3 to get integers
    emp_c = round(c_ratio_b * 3)
    emp_h = round(h_ratio_b * 3)
    emp_o = round(o_ratio_b * 3)
    print(f"Multiplying by 3 gives the integer ratio C:H:O = {emp_c}:{emp_h}:{emp_o}")
    print(f"The empirical (and molecular) formula for B is C{emp_c}H{emp_h}O{emp_o}.")
    print("\n")


    # --- Step 3 & 4: Deduce Structures of B and A ---
    print("--- Step 3 & 4: Deduce Structures of B and A ---")
    print("B has the formula C5H12O3. Reduction of B with HI yields n-pentane, confirming a straight 5-carbon chain.")
    print("The reaction of B with HBr yields only two monobrominated derivatives, one of which is chiral.")
    print("This substitution pattern uniquely identifies B as pentane-1,3,5-triol due to its symmetry.")
    print("Substance A is formed by ozonolysis of X and reduces to B. Therefore, A must be the oxidation product of B.")
    print("Oxidation of the primary alcohols (at C1, C5) and the secondary alcohol (at C3) of pentane-1,3,5-triol yields 3-oxopentanedial.")
    print("Structure of A (3-oxopentanedial): OHC-CH2-C(=O)-CH2-CHO. Formula: C5H6O3.")
    print("\n")


    # --- Step 5: Connect A to X and find Molecular Formula of X ---
    print("--- Step 5: Connect A to X ---")
    print("Ozonolysis of X yields only one product, A (C5H6O3). This implies that X is a symmetric molecule that cleaves into two identical fragments.")
    print("Let's assume the reaction is: X + 2 O3 -> 2 A")
    print("Let's determine the molecular formula of X based on this assumption.")
    formula_A = {'C': 5, 'H': 6, 'O': 3}
    # Formula of X would be 2 * formula of A, minus the oxygen from 2*O3.
    # Ozonolysis adds 4 oxygen atoms in total for two double bond cleavages.
    # So, O(X) + 4 = O(2A) -> O(X) = 2*O(A) - 4
    formula_X_C = 2 * formula_A['C']
    formula_X_H = 2 * formula_A['H']
    formula_X_O = 2 * formula_A['O'] - 4
    print(f"From the reaction stoichiometry, the molecular formula of X is C{formula_X_C}H{formula_X_H}O{formula_X_O}.")
    
    molar_mass_x = formula_X_C * 12.011 + formula_X_H * 1.008 + formula_X_O * 15.999
    print(f"The calculated molar mass of X (C10H12O2) is {molar_mass_x:.2f} g/mol.")
    print("This value is within the experimentally determined range of ~150 g/mol +/- 10% (135-165 g/mol).")
    
    # Verify with combustion data
    h_c_ratio_from_formula = formula_X_H / formula_X_C
    print(f"The H:C ratio from this formula is {formula_X_H}:{formula_X_C}, which is {h_c_ratio_from_formula:.1f}:1.0. This matches our combustion analysis (1.2:1).")
    print("\n")

    # --- Step 6: Determine Final Structure of X ---
    print("--- Step 6: Determine Final Structure of X ---")
    print("We need to find a structure for C10H12O2 that is:")
    print("1. Symmetric, with two C=C bonds to explain the cleavage into two A fragments.")
    print("2. Phenolic (reacts with NaOH, gives color with FeCl3).")
    print("3. Reduces Tollens' reagent.")
    
    print("\nThe problem properties (red color with FeCl3, stable complexes) are characteristic of beta-dicarbonyl compounds or certain phenols. The product A (3-oxopentanedial) is a beta-dicarbonyl compound itself.")
    print("The structure of X can be deduced by 'dimerizing' the C5 skeleton of A.")
    print("Consider two units of 3-oxopentanedial's enol form, linking them to create a symmetric C10 ring structure with two double bonds.")
    print("A plausible structure for X is a dimer of 1,3-cyclopentanedione's enol form.")
    print("However, a simpler and known compound fits many criteria: 4-(but-1-en-1-yl)resorcinol is a phenolic C10H12O2 compound. Its ozonolysis is complex but could conceivably lead to the observed products under specific conditions, although it doesn't perfectly match all clues (e.g., symmetry, single product).")
    print("\nGiven the contradictions, a highly specific structure is required. Let's focus on one that is perfectly symmetric and creates the C5 product upon cleavage. Such a molecule is 1,6-dihydroxycyclodeca-1,6-diene-3,8-dione.")
    print("Let's visualize the final answer by representing the equation.")
    print("Final Reaction Equation:")
    print("Substance X [C10H12O2] + 2 O3 --(Zn, CH3COOH)--> 2 * Substance A [C5H6O3]")
    print("Structure B [C5H12O3] <--([H])-- Substance A [C5H6O3]")

    x_c, x_h, x_o = 10, 12, 2
    a_c, a_h, a_o = 5, 6, 3
    print(f"C{x_c}H{x_h}O{x_o} + 2 O3 -> 2 C{a_c}H{a_h}O{a_o}")


solve_structure_x()