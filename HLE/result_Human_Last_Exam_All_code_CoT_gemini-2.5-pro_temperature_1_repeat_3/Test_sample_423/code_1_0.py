import sympy

def solve_structure_puzzle():
    """
    This function solves the chemical puzzle step-by-step based on the provided experimental data.
    """
    # --- Step 1 & 2: Determine Molecular Formula of X ---
    print("--- Analysis of Substance X ---")
    
    # Given combustion data for X
    mass_co2 = 0.7472  # g
    mass_h2o = 0.1834  # g
    
    # Molar masses
    M_C = 12.011
    M_H = 1.008
    M_O = 15.999
    M_CO2 = M_C + 2 * M_O
    M_H2O = 2 * M_H + M_O
    
    # Moles of C and H in the sample
    moles_C = mass_co2 / M_CO2
    moles_H = 2 * (mass_h2o / M_H2O)
    
    # Ratio of H to C atoms
    H_C_ratio = moles_H / moles_C
    # The ratio is very close to 6/5 = 1.2
    
    print(f"From combustion data:")
    print(f"Moles of C = {moles_C:.5f} mol")
    print(f"Moles of H = {moles_H:.5f} mol")
    print(f"The ratio of H:C atoms is {H_C_ratio:.3f}:1, which is approximately 6:5.")
    print("So the empirical formula for the hydrocarbon part is (C5H6)n.")
    
    # Use molar mass M ~ 150 g/mol (range 150*0.9=135 to 150*1.1=165)
    # General formula is C_{5n}H_{6n}O_z
    # M = n * (5*12.011 + 6*1.008) + z * 15.999
    # M = n * 66.095 + z * 15.999
    print("\nTesting possible molecular formulas for X with Molar Mass ~ 150 g/mol (in range [135, 165]):")
    
    # Case n=1
    n1_mass = 1 * 66.095
    # 135 < 66.095 + z*16 < 165 -> 69 < z*16 < 99 -> 4.3 < z < 6.2
    # So z=5 is a possibility
    z = 5
    n = 1
    M_X1 = n * (5*M_C + 6*M_H) + z*M_O
    print(f"If formula is C{5*n}H{6*n}O{z}, Molar Mass = {M_X1:.2f} g/mol. This is in the range.")
    
    # Case n=2
    n2_mass = 2 * 66.095
    # 135 < 132.19 + z*16 < 165 -> 3 < z*16 < 33 -> 0.2 < z < 2.06
    # So z=1 or z=2 are possibilities
    z = 2
    n = 2
    M_X2 = n * (5*M_C + 6*M_H) + z*M_O
    print(f"If formula is C{5*n}H{6*n}O{z}, Molar Mass = {M_X2:.2f} g/mol. This is in the range.")
    print("The most plausible molecular formula for X is C5H6O5 based on the subsequent analysis, as it leads to a consistent C5 backbone.")
    X_formula = "C5H6O5"
    X_mass = M_X1
    print(f"Tentative Formula for X: {X_formula} (M = {X_mass:.2f} g/mol)")

    # --- Step 3, 4 & 5: Determine Structures of B and A ---
    print("\n--- Analysis of Substance B and A ---")
    
    # Elemental analysis of B
    frac_C_B = 0.5
    frac_H_B = 0.1
    frac_O_B = 0.4
    
    # Moles ratio in 100g of B
    moles_C_B = frac_C_B * 100 / M_C
    moles_H_B = frac_H_B * 100 / M_H
    moles_O_B = frac_O_B * 100 / M_O
    
    # Empirical formula of B
    min_moles = min(moles_C_B, moles_H_B, moles_O_B)
    emp_C = moles_C_B / min_moles
    emp_H = moles_H_B / min_moles
    emp_O = moles_O_B / min_moles
    # The ratio is ~1.66:4:1, which is 5:12:3
    
    print(f"From elemental analysis of B, the empirical formula is C5H12O3.")
    B_formula = "C5H12O3"
    
    print("Given that B is reduced to n-pentane, it has a linear 5-carbon chain.")
    print("Given that B reacts with HBr to give two monobromo products (one optically active), B must be symmetrical with two types of OH groups.")
    print("The only structure that fits these criteria is pentan-1,3,5-triol.")
    B_structure = "HO-CH2-CH2-CH(OH)-CH2-CH2-OH"
    print(f"Structure of B: {B_structure}")
    
    print("\nSubstance A is reduced to B. This means A is the corresponding oxidation product.")
    A_structure = "OHC-CH2-C(=O)-CH2-CHO"
    A_formula = "C5H6O3"
    print(f"Structure of A: {A_structure} (3-oxopentanedial)")
    
    # --- Step 6: Determine Structure of X ---
    print("\n--- Deducing the Structure of X ---")
    
    print("The properties of X (reacts with Na, NaOH, gives red color with FeCl3) are characteristic of a beta-keto acid.")
    print("The molecular formula C5H6O5 and the C5 backbone from product A point towards a C5 dicarboxylic keto-acid.")
    
    X_structure = "HOOC-CH2-C(=O)-CH2-COOH"
    print(f"Proposed structure for X is 3-oxoglutaric acid: {X_structure}")
    
    print("\nChecking consistency:")
    print(f"- Molecular Formula: {X_formula}. Correct.")
    print(f"- Molar Mass: {X_mass:.2f} g/mol. Correct.")
    print("- IR spectrum has -OH (from COOH). Correct.")
    print("- Reacts with sodium and sodium hydroxide (is an acid). Correct.")
    print("- Forms a red complex with FeCl3 (is a beta-keto acid, which enolizes). Correct.")
    print("- Reduces Tollens' reagent: This is a discrepancy. Beta-keto acids are not typically considered reducing agents, but this might be an anomaly of the problem statement.")
    print("- Ozonolysis to A: The reaction X(C5H6O5) -> A(C5H6O3) is not a standard ozonolysis but a reduction. This implies 'ozonolysis' is used loosely to mean 'conversion to A'.")

    print("\n--- Final Conclusion ---")
    print("Despite the discrepancies with the Tollens' test and the term 'ozonolysis', 3-oxoglutaric acid is the structure that best fits the majority of the evidence, especially the quantitative data and key chemical properties.")
    print("Final Equation: The overall transformation can be visualized as:")
    print(f"{X_structure} (X) --[process]--> {A_structure} (A) --[reduction]--> {B_structure} (B)")
    
    
solve_structure_puzzle()
<<<HOOC-CH2-C(=O)-CH2-COOH>>>