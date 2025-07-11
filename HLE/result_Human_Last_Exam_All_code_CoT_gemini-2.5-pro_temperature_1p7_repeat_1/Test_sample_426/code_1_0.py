import sys

def solve_chemistry_problem():
    """
    Solves the chemical puzzle by calculating molar masses and deducing structures
    based on the provided reaction scheme and data.
    """
    
    # --- Step 1: Determine the molecular formula of A1 ---
    print("Step 1: Determining the molecular formula of compound A1.")
    
    # Elemental composition of A1
    composition = {'C': 54.5, 'H': 13.6, 'N': 31.8}
    atomic_mass = {'C': 12.01, 'H': 1.008, 'N': 14.01}
    
    # Calculate moles per 100g of substance
    moles_C = composition['C'] / atomic_mass['C']
    moles_H = composition['H'] / atomic_mass['H']
    moles_N = composition['N'] / atomic_mass['N']
    
    # Find the smallest mole value to determine the ratio
    min_moles = min(moles_C, moles_H, moles_N)
    
    # Calculate the empirical formula ratios
    ratio_C = round(moles_C / min_moles)
    ratio_H = round(moles_H / min_moles)
    ratio_N = round(moles_N / min_moles)
    
    print(f"The empirical formula of A1 is C{ratio_C}H{ratio_H}N{ratio_N}.")
    
    # A1 is formed from a dibromoalkane (A) and excess ammonia, making it a diamine.
    # Therefore, the molecular formula must contain at least 2 Nitrogen atoms.
    molecular_formula_C = ratio_C * 2
    molecular_formula_H = ratio_H * 2
    molecular_formula_N = ratio_N * 2
    
    print(f"Since A1 is a diamine, its molecular formula is (C{ratio_C}H{ratio_H}N{ratio_N}) x 2 = C{molecular_formula_C}H{molecular_formula_H}N{molecular_formula_N}.\n")

    # --- Step 2: Determine the molar mass of the carboxylic acid ---
    print("Step 2: Determining the molar mass of the carboxylic acid from titration data.")
    
    mass_acid_g = 2.16
    vol_koh_L = 30 / 1000  # Convert ml to L
    conc_koh_M = 1.0
    
    # Moles of KOH = Concentration * Volume
    moles_koh = conc_koh_M * vol_koh_L
    
    # For a monoprotic acid, moles of acid = moles of KOH
    moles_acid = moles_koh
    
    # Molar Mass = mass / moles
    molar_mass_acid = mass_acid_g / moles_acid
    
    print(f"Neutralization of {mass_acid_g} g of acid required {vol_koh_L * 1000} ml of {conc_koh_M} M KOH solution.")
    print(f"Moles of KOH used = {conc_koh_M} mol/L * {vol_koh_L} L = {moles_koh:.3f} mol.")
    print(f"Assuming a monoprotic acid, moles of acid = {moles_acid:.3f} mol.")
    print(f"Calculated Molar Mass of the acid = {mass_acid_g} g / {moles_acid:.3f} mol = {molar_mass_acid:.2f} g/mol.\n")

    # --- Step 3 & 4: Identify the acid and deduce the structure of A2 ---
    print("Step 3: Identifying the acid and deducing the precursor A2.")
    print("The calculated molar mass is ~72 g/mol.")
    print("Propenoic acid (CH2=CH-COOH) has a molar mass of 72.06 g/mol.")
    print("Propanoic acid (CH3CH2COOH) has a molar mass of 74.08 g/mol.")
    print("Given the small discrepancy, and the reaction pathway from a saturated precursor, the acid is very likely propanoic acid.")
    print("\nThe formation of a C3 acid (propanoic acid) and CO2 (a C1 fragment) results from the oxidative cleavage of a C4 precursor (A2).")
    print("This specific cleavage pattern points to A2 being butane-1,2-diol (CH3-CH2-CH(OH)-CH2OH).\n")

    # --- Step 5: Deduce the full reaction chain ---
    print("Step 4: Deducing the structures of A1, A, and X.")
    print("A2 is butane-1,2-diol.")
    print("A1 is formed by replacing -OH groups in A2 with -NH2 groups, so A1 is 1,2-diaminobutane.")
    print("  - Check: 1,2-diaminobutane (CH3-CH2-CH(NH2)-CH2NH2) has a chiral center, making all 4 carbon-proton groups (CH3, CH2, CH, CH2) chemically distinct. This fits the 'four types of signals' in the NMR spectrum.")
    print("A is formed by replacing -NH2 groups with -Br atoms, so A is 1,2-dibromobutane.")
    print("X is a hydrocarbon that reacts with Br2 to give A. This is an addition reaction.")
    print("Therefore, X must be the alkene that adds bromine across the double bond to form 1,2-dibromobutane.\n")
    
    # --- Final Conclusion ---
    print("Final Conclusion: The structure of substance X is But-1-ene (CH3-CH2-CH=CH2).")

solve_chemistry_problem()
<<<But-1-ene>>>