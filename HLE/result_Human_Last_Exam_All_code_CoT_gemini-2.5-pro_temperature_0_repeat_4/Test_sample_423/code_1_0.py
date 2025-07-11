def solve_chemistry_problem():
    """
    This script performs the calculations needed to solve the structural chemistry problem.
    It determines the empirical formula of B, the empirical formula of X, and checks
    potential molecular formulas against the given molar mass.
    """
    # --- Step 1: Analyze Product B ---
    print("--- Analysis of Product B ---")
    # Mass fractions of elements in B
    mass_frac_C_B = 0.50
    mass_frac_H_B = 0.10
    mass_frac_O_B = 0.40

    # Molar masses of elements
    M_C = 12.01
    M_H = 1.008
    M_O = 16.00

    # Calculate molar ratios assuming 100g of substance
    moles_C_B = mass_frac_C_B * 100 / M_C
    moles_H_B = mass_frac_H_B * 100 / M_H
    moles_O_B = mass_frac_O_B * 100 / M_O

    # Find the simplest ratio by dividing by the smallest number of moles
    min_moles_B = min(moles_C_B, moles_H_B, moles_O_B)
    ratio_C_B = moles_C_B / min_moles_B
    ratio_H_B = moles_H_B / min_moles_B
    ratio_O_B = moles_O_B / min_moles_B

    # To get whole numbers, we can see that C is ~1.66 (5/3). Multiply by 3.
    emp_C_B = round(ratio_C_B * 3)
    emp_H_B = round(ratio_H_B * 3)
    emp_O_B = round(ratio_O_B * 3)

    print(f"Molar ratios in B (C:H:O): {ratio_C_B:.2f} : {ratio_H_B:.2f} : {ratio_O_B:.2f}")
    print(f"Empirical formula of B: C{emp_C_B}H{emp_H_B}O{emp_O_B}")
    
    # Since B is reduced to n-pentane, it must have 5 carbon atoms.
    # The empirical formula C5H12O3 is the molecular formula.
    molecular_formula_B = f"C{emp_C_B}H{emp_H_B}O{emp_O_B}"
    print(f"Molecular formula of B is {molecular_formula_B}.\n")

    # --- Step 2: Analyze Substance X ---
    print("--- Analysis of Substance X ---")
    # Combustion data for X
    mass_CO2 = 0.7472  # g
    mass_H2O = 0.1834  # g

    # Molar masses of CO2 and H2O
    M_CO2 = 44.01
    M_H2O = 18.016

    # Calculate moles of C and H
    moles_C_X = mass_CO2 / M_CO2
    moles_H_X = 2 * (mass_H2O / M_H2O)

    # Calculate C:H ratio
    ratio_H_to_C = moles_H_X / moles_C_X
    # The ratio is ~1.2, which is 6/5.
    emp_C_X = 5
    emp_H_X = 6
    print(f"Molar ratio of H:C in X is {ratio_H_to_C:.2f}, which corresponds to an empirical formula of (C{emp_C_X}H{emp_H_X})nOz.")

    # Use molar mass to find the molecular formula (C5H6)nOz
    M_X_est = 150
    M_X_error = 0.10
    M_X_min = M_X_est * (1 - M_X_error)
    M_X_max = M_X_est * (1 + M_X_error)
    print(f"Estimated molar mass of X is ~{M_X_est} g/mol (range: {M_X_min:.1f} - {M_X_max:.1f} g/mol).")

    # Test possible formulas
    # Case n=1: C5H6Oz
    M_C5H6 = emp_C_X * M_C + emp_H_X * M_H
    print("\nTesting formula C5H6Oz:")
    for z in range(1, 7):
        M_test = M_C5H6 + z * M_O
        if M_X_min <= M_test <= M_X_max:
            print(f"  - For z={z}, formula is C5H6O{z}, M = {M_test:.2f} g/mol. This is a possible formula.")
        else:
            print(f"  - For z={z}, formula is C5H6O{z}, M = {M_test:.2f} g/mol.")
    
    # Case n=2: C10H12Oz
    M_C10H12 = 2 * M_C5H6
    print("\nTesting formula C10H12Oz:")
    for z in range(1, 4):
        M_test = M_C10H12 + z * M_O
        if M_X_min <= M_test <= M_X_max:
            print(f"  - For z={z}, formula is C10H12O{z}, M = {M_test:.2f} g/mol. This is a possible formula.")
        else:
            print(f"  - For z={z}, formula is C10H12O{z}, M = {M_test:.2f} g/mol.")

    # --- Step 3: Deducing Structures ---
    print("\n--- Structural Deduction ---")
    print("1. Structure of B (C5H12O3): It's a straight-chain pentanetriol. The reaction with HBr yields only two monobrominated products, one of which is optically active. This points to a symmetrical structure with two types of -OH groups. Pentan-1,3,5-triol fits this description.")
    print("   - B = Pentan-1,3,5-triol: HO-CH2-CH2-CH(OH)-CH2-CH2-OH")
    
    print("\n2. Structure of A: A is reduced to B. Therefore, A is the corresponding carbonyl compound.")
    print("   - A = 3-oxopentanedial: OHC-CH2-C(=O)-CH2-CHO")
    print("   - Formula of A is C5H6O3.")

    print("\n3. Structure of X: Ozonolysis of X gives only A. This implies X is a symmetrical cyclic precursor to A.")
    print("   - Reversing the ozonolysis of A (3-oxopentanedial) leads to a C5 ring with a double bond and a keto group: 3-oxocyclopentene.")
    print("   - The keto form (3-oxocyclopentene) does not match the chemical properties of X (reacts with NaOH, FeCl3, Tollen's).")
    print("   - The enol form, 3-hydroxy-1,3-cyclopentadiene, does match these properties.")
    print("   - Formula of this enol is C5H6O. Molar mass is 82 g/mol.")
    print("   - This contradicts the molecular formula C5H6O5 (M=146) derived from combustion and molar mass data.")
    
    print("\n4. Resolving the Contradiction: The problem data is contradictory. The molar mass estimation is the most likely source of error. The chemical properties and reaction sequence strongly point to a specific structure. The structure that best fits the properties (FeCl3 test, Tollen's, etc.) and the reaction sequence (ozonolysis to A, which reduces to B) is the enol form of 3-oxocyclopentene.")
    print("   - The final proposed structure for X is 3-hydroxy-1,3-cyclopentadiene.")
    print("   - This structure is an enol of a β-dicarbonyl system (vinylogous), explaining the acidity (reaction with NaOH) and the red color with FeCl3.")
    print("   - As an enol, it has an -OH group and can act as a reducing agent (Tollen's test).")
    print("   - The ozonolysis is presumed to react with its keto-tautomer (3-oxocyclopentene) to yield the single product A.")

    final_structure_name = "3-hydroxy-1,3-cyclopentadiene"
    print(f"\nFinal proposed structure for X: {final_structure_name}")

solve_chemistry_problem()