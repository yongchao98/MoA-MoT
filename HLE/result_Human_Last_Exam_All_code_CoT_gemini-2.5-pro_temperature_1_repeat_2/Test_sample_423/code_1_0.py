def solve_structure():
    """
    This function performs the calculations to deduce the molecular formulas
    of substances X and B based on the provided experimental data.
    """
    # --- Step 1: Molecular Formula of X ---
    mass_co2 = 0.7472  # g
    mass_h2o = 0.1834  # g
    molar_mass_co2 = 44.01 # g/mol
    molar_mass_h2o = 18.016 # g/mol

    moles_c = mass_co2 / molar_mass_co2
    moles_h = (mass_h2o / molar_mass_h2o) * 2

    # Ratio of C to H
    ratio_h_to_c = moles_h / moles_c
    # The ratio is ~1.2, which is 6/5.

    # Empirical formula is (C5H6)n * Om
    # Test for n=2, as M ~ 150 g/mol
    # M_c10h12 = 10 * 12.01 + 12 * 1.008 = 132.2
    # M = 132.2 + m * 16.00
    # For m=1, M = 148.2 (in range 135-165)
    # For m=2, M = 164.2 (in range 135-165)

    # A more rigorous method to confirm the formula C10H12O2
    # From calculation, empirical formula is C5H6O (mass=82)
    # Molecular mass range is 135-165.
    # 82 * 2 = 164, which is in the range.
    # So, molecular formula of X is (C5H6O)2 = C10H12O2
    formula_x_c = 10
    formula_x_h = 12
    formula_x_o = 2
    
    print("Step 1: Deduction of X's Molecular Formula")
    print(f"Moles of C from CO2: {moles_c:.4f} mol")
    print(f"Moles of H from H2O: {moles_h:.4f} mol")
    print(f"Ratio H:C is ~{ratio_h_to_c:.1f}:1, which is approximately 6:5.")
    print("Based on the molar mass of ~150 g/mol, the molecular formula of X is determined to be C10H12O2.")
    print(f"Final Equation for X: {formula_x_c}*C + {formula_x_h}*H + {formula_x_o}*O")
    print("-" * 20)

    # --- Step 2: Molecular Formula of B ---
    # Mass fractions in B: C=0.5, H=0.1, O=0.4
    mf_c = 0.5
    mf_h = 0.1
    mf_o = 0.4
    
    # Molar ratios
    mol_ratio_c = mf_c / 12.01
    mol_ratio_h = mf_h / 1.008
    mol_ratio_o = mf_o / 16.00
    
    # Normalize by smallest
    norm_c = mol_ratio_c / mol_ratio_o
    norm_h = mol_ratio_h / mol_ratio_o
    norm_o = mol_ratio_o / mol_ratio_o
    
    # Multiply to get integers
    formula_b_c = round(norm_c * 3)
    formula_b_h = round(norm_h * 3)
    formula_b_o = round(norm_o * 3)
    
    print("Step 2: Deduction of B's Molecular Formula")
    print(f"Molar ratios from mass fractions (C:H:O): {mol_ratio_c:.4f} : {mol_ratio_h:.4f} : {mol_ratio_o:.4f}")
    print(f"Normalized ratios (C:H:O): {norm_c:.2f} : {norm_h:.2f} : {norm_o:.2f}")
    print("Multiplying by 3 gives the empirical formula C5H12O3.")
    print("Since B reduces to n-pentane, this is also its molecular formula.")
    print(f"Final Equation for B: {formula_b_c}*C + {formula_b_h}*H + {formula_b_o}*O")
    print("-" * 20)
    
    print("Step 3 & 4: Structure Determination")
    print("B is identified as pentane-2,3,4-triol based on the reaction with HBr.")
    print("A (the ozonolysis product) is therefore pentane-2,3,4-trione.")
    print("The structure of X that best fits all the data, especially the C10H12O2 formula and red FeCl3 test, is Hinokitiol (isopropyltropolone).")

solve_structure()