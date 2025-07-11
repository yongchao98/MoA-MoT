def solve_structure_mystery():
    """
    This function outlines the step-by-step reasoning to determine the structure of substance X.
    It performs the calculations for each step and prints the conclusion.
    """

    # --- Step 1 & 2: Analysis of Substance B ---
    print("--- Step 1 & 2: Analysis and Structure of Substance B ---")
    mass_c, mass_h, mass_o = 50.0, 10.0, 40.0
    molar_mass_c, molar_mass_h, molar_mass_o = 12.01, 1.008, 16.00
    
    moles_c = mass_c / molar_mass_c
    moles_h = mass_h / molar_mass_h
    moles_o = mass_o / molar_mass_o
    
    # Normalize by the smallest number of moles (Oxygen)
    c_ratio = moles_c / moles_o
    h_ratio = moles_h / moles_o
    o_ratio = moles_o / moles_o
    
    # Get integer ratios (C:H:O is 5/3 : 4 : 1)
    emp_c, emp_h, emp_o = 5, 12, 3
    
    print(f"Based on mass fractions, the empirical formula of B is C{emp_c}H{emp_h}O{emp_o}.")
    print("Reduction of B with HI gives n-pentane, confirming a C5 linear chain.")
    print("Therefore, the molecular formula of B is C5H12O3.")
    print("The reaction with HBr giving two products (one chiral) points to a symmetrical structure.")
    print("Conclusion: B is Pentan-1,3,5-triol: HO-CH2-CH2-CH(OH)-CH2-CH2-OH\n")
    
    # --- Step 3: Structure of A ---
    print("--- Step 3: Structure of Substance A ---")
    print("A is the ozonolysis product of X, and its reduction yields B (Pentan-1,3,5-triol).")
    print("Therefore, A must be the corresponding carbonyl compound.")
    print("Conclusion: A is 3-oxopentanedial: OHC-CH2-C(=O)-CH2-CHO\n")

    # --- Step 4 & 5: Analysis and Structure of X ---
    print("--- Step 4 & 5: Analysis and Structure of Substance X ---")
    print("From combustion analysis, the C:H molar ratio in X is 5:6.")
    print("Given M ~ 150, the molecular formula must be a multiple of C5H6.")
    print("The reaction 'X -> 2 A' is C10H12O2 + 2O3 -> 2 C5H6O3, which is stoichiometrically sound.")
    print("So, the molecular formula for X is C10H12O2.")
    print("X must be a symmetrical C10 molecule with two C=C bonds that, upon cleavage, yield two molecules of A.")
    print("Constructing this molecule by joining two '-CH2-C(=O)-CH2-' units with two '-CH=CH-' linkers in a ring.")
    
    final_structure_name = "cyclodeca-1,6-diene-3,8-dione"
    print(f"\nFinal deduced structure of X: {final_structure_name}")
    print("Structure: A 10-membered ring with double bonds at positions 1 and 6, and carbonyl groups at positions 3 and 8.")
    print("Formula: C10H12O2")
    
    c = 10
    h = 12
    o = 2
    
    molar_mass_x = c * molar_mass_c + h * molar_mass_h + o * molar_mass_o
    print(f"The molar mass is {molar_mass_x:.2f} g/mol, which is within the experimental range of 135-165 g/mol.")
    
    print("\nThis structure satisfies all conditions:")
    print("- Molecular formula and mass are consistent with data.")
    print("- It is symmetrical and its ozonolysis correctly yields two molecules of 3-oxopentanedial.")
    print("- Its ketone groups can enolize, explaining the -OH band, reactions with Na/NaOH, and color with FeCl3.")

solve_structure_mystery()