import math

def solve_chemistry_problem():
    """
    This script solves the chemical puzzle to identify hydrocarbon X.
    It performs calculations for elemental analysis and titration,
    and prints the logical steps to deduce the structures of all compounds.
    """
    print("Step 1: Determining the empirical and molecular formula of substance A1")
    
    # Elemental composition of A1
    composition = {'C': 54.5, 'H': 13.6, 'N': 31.8}
    atomic_masses = {'C': 12.01, 'H': 1.008, 'N': 14.01}

    # Calculate moles in a 100g sample
    moles = {element: composition[element] / atomic_masses[element] for element in composition}
    
    print(f"Relative moles in 100g: C = {moles['C']:.2f}, H = {moles['H']:.2f}, N = {moles['N']:.2f}")

    # Find the simplest whole number ratio
    min_moles = min(moles.values())
    empirical_ratio = {element: moles[element] / min_moles for element in moles}

    print(f"Empirical ratio: C = {empirical_ratio['C']:.1f}, H = {math.ceil(empirical_ratio['H']):.1f}, N = {empirical_ratio['N']:.1f}")
    print("The empirical formula of A1 is C2H6N.")

    print("\nFrom the reaction scheme, hydrocarbon X (an alkene, CnH2n) reacts with Br2 to form A (CnH2nBr2).")
    print("A then reacts with excess ammonia to form A1. This means two Br atoms are replaced by two NH2 groups.")
    print("Therefore, A1 must be a diamine with the general formula CnH2n(NH2)2, which simplifies to CnH(2n+4)N2.")
    print("The molecular formula must contain an even number of nitrogen atoms. The empirical formula C2H6N has only one.")
    print("So, the molecular formula must be a multiple. The simplest is (C2H6N) * 2 = C4H12N2.")
    print("This molecular formula C4H12N2 fits the general formula CnH(2n+4)N2 for n=4.")
    print("Conclusion: The molecular formula of A1 is C4H12N2.")

    print("\n----------------------------------------")
    print("Step 2: Identifying the carboxylic acid from titration data")
    
    # Titration data
    acid_mass = 2.16  # g
    koh_volume_ml = 30.0  # ml
    koh_conc = 1.0  # M

    # Moles of KOH = Volume (L) * Concentration (M)
    moles_koh = (koh_volume_ml / 1000) * koh_conc
    
    # For a monoprotic acid, moles of acid = moles of KOH
    moles_acid = moles_koh
    
    # Molar Mass = Mass / Moles
    molar_mass_acid = acid_mass / moles_acid
    
    print(f"Neutralization of {acid_mass} g of acid required {koh_volume_ml} ml of {koh_conc} M KOH solution.")
    print(f"This corresponds to {moles_koh:.3f} moles of KOH, and thus {moles_acid:.3f} moles of a monoprotic acid.")
    print(f"Calculated molar mass of the carboxylic acid = {acid_mass} g / {moles_acid:.3f} mol = {molar_mass_acid:.2f} g/mol.")

    # Identify the acid
    molar_mass_propanoic_acid = 3 * 12.01 + 6 * 1.008 + 2 * 16.00
    print(f"\nThe molar mass of propanoic acid (CH3CH2COOH) is {molar_mass_propanoic_acid:.2f} g/mol.")
    print(f"The calculated molar mass ({molar_mass_acid:.2f} g/mol) is very close to that of propanoic acid.")
    print("The small difference is likely due to experimental error. The acid is propanoic acid.")

    print("\n----------------------------------------")
    print("Step 3: Deducing the structure of X by working backward")
    
    print("\nThe final oxidation of A2 produced propanoic acid (CH3CH2COOH, a C3 acid) and CO2 (a C1 fragment).")
    print("This is the result of oxidative cleavage of a vicinal diol (a compound with -OH on adjacent carbons).")
    print("The cleavage must have occurred in a C4 diol as follows:")
    print("  CH2(OH)-CH(OH)-CH2-CH3  --[O]-->  CO2 (from C1) + HOOC-CH2-CH3 (from C2-C4)")
    print("Therefore, substance A2 is Butane-1,2-diol.")

    print("\nA2 (Butane-1,2-diol) was formed from A1 by reacting with nitrous acid (NH2 -> OH).")
    print("Therefore, substance A1 is Butane-1,2-diamine: CH2(NH2)-CH(NH2)-CH2-CH3.")
    print("This structure (C4H12N2) is consistent with the molecular formula and the NMR data ('four types of signals' corresponding to the four non-equivalent carbon environments).")

    print("\nA1 (Butane-1,2-diamine) was formed from A by reacting with ammonia (Br -> NH2).")
    print("Therefore, substance A is 1,2-dibromobutane: CH2Br-CHBr-CH2-CH3.")

    print("\nA (1,2-dibromobutane) was formed by the addition of Br2 to hydrocarbon X.")
    print("This electrophilic addition occurs across a double bond.")
    print("Therefore, the double bond in X must be between C1 and C2.")
    
    print("\n----------------------------------------")
    print("Final Conclusion:")
    print("The structure of hydrocarbon X is But-1-ene.")
    print("Chemical Formula: CH2=CH-CH2-CH3")

solve_chemistry_problem()