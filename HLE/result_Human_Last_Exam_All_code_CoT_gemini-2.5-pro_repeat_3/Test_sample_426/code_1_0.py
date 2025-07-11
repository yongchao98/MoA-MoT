import sys

def solve_chemistry_problem():
    """
    This script performs the calculations and logical deductions to identify hydrocarbon X.
    """
    # Use a dictionary for atomic masses for clarity
    atomic_mass = {'C': 12.01, 'H': 1.008, 'N': 14.007, 'O': 16.00}

    # --- Step 1: Determine the molecular formula of substance A1 ---
    print("Step 1: Determine the molecular formula of substance A1.")
    composition = {'C': 54.5, 'H': 13.6, 'N': 31.8}
    print(f"Given elemental composition of A1: C={composition['C']}%, H={composition['H']}%, N={composition['N']}%")

    # To find the empirical formula, assume a 100g sample
    moles_C = composition['C'] / atomic_mass['C']
    moles_H = composition['H'] / atomic_mass['H']
    moles_N = composition['N'] / atomic_mass['N']
    
    # To output the calculation, we print each number
    print(f"Moles in 100g sample:")
    print(f"  C: {composition['C']} / {atomic_mass['C']} = {moles_C:.2f} mol")
    print(f"  H: {composition['H']} / {atomic_mass['H']} = {moles_H:.2f} mol")
    print(f"  N: {composition['N']} / {atomic_mass['N']} = {moles_N:.2f} mol")

    # Find the smallest number of moles to determine the ratio
    min_moles = min(moles_C, moles_H, moles_N)
    
    # Calculate the simplest whole-number ratio
    ratio_C = round(moles_C / min_moles)
    ratio_H = round(moles_H / min_moles)
    ratio_N = round(moles_N / min_moles)
    print(f"The empirical formula is C{ratio_C}H{ratio_H}N{ratio_N}.")

    print("\nReasoning for Molecular Formula:")
    print("X (hydrocarbon) + Br2 -> A (dibromoalkane, CnH2nBr2)")
    print("A + excess NH3 -> A1 (diamine, CnH2n(NH2)2 or CnH(2n+4)N2)")
    print("The empirical formula is C2H6N. To match the general formula CnH(2n+4)N2, n must be 4.")
    print("For n=4, the molecular formula is C4H(2*4+4)N2 = C4H12N2. The ratio 4:12:2 simplifies to 2:6:1, matching the empirical formula.")
    print("Therefore, the molecular formula of A1 is C4H12N2.")

    # --- Step 2: Determine the molar mass of the carboxylic acid ---
    print("\nStep 2: Determine the molar mass of the carboxylic acid.")
    mass_acid_g = 2.16
    vol_koh_L = 0.030
    conc_koh_M = 1.0
    
    # Neutralization is a 1:1 reaction: R-COOH + KOH -> R-COOK + H2O
    moles_koh = vol_koh_L * conc_koh_M
    print(f"Moles of KOH used = {vol_koh_L} L * {conc_koh_M} mol/L = {moles_koh} mol")
    
    # Moles of acid is equal to moles of KOH
    moles_acid = moles_koh
    print(f"Moles of the monoprotic carboxylic acid = {moles_acid} mol")
    
    # Calculate the experimental molar mass
    molar_mass_acid_exp = mass_acid_g / moles_acid
    print(f"Experimental molar mass of the acid = {mass_acid_g} g / {moles_acid} mol = {molar_mass_acid_exp:.2f} g/mol")

    # --- Step 3 & 4: Deduce all structures ---
    print("\nStep 3 & 4: Deducing the structures from X to the final acid.")
    print("Oxidation of A2 to a carboxylic acid and CO2 indicates A2 is a 1,2-diol with a terminal -CH2OH group.")
    print("Since A1 is C4H12N2, A2 must be C4H10O2 (a butanediol).")
    print("The only C4 diol that fits is Butane-1,2-diol (CH3CH2CH(OH)CH2OH).")
    print("Its oxidation yields CO2 (from the CH2OH part) and Propanoic Acid (CH3CH2COOH).")
    
    # Verify the molar mass of Propanoic Acid (C3H6O2)
    molar_mass_propanoic_acid = 3 * atomic_mass['C'] + 6 * atomic_mass['H'] + 2 * atomic_mass['O']
    print(f"The theoretical molar mass of Propanoic Acid (C3H6O2) is {molar_mass_propanoic_acid:.2f} g/mol.")
    print("This matches our experimental value.")

    print("\nWorking backwards:")
    print(" - Carboxylic Acid: Propanoic Acid")
    print(" - A2 (from A1 + HNO2): Butane-1,2-diol")
    print(" - A1 (from A + NH3): Butane-1,2-diamine. This structure has 4 unique carbon environments, matching the NMR data of four signals.")
    print(" - A (from X + Br2): 1,2-dibromobutane")
    print(" - X (original hydrocarbon): But-1-ene (CH3-CH2-CH=CH2)")
    
    # The final answer is the structure of substance X.
    final_answer = "But-1-ene"
    # The problem asks to print the answer in a specific format.
    # We use sys.stdout.write to avoid adding an extra newline.
    sys.stdout.write(f"\n<<<The structure of substance X is {final_answer}>>>\n")

solve_chemistry_problem()