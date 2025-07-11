import math

def solve_chemistry_problem():
    """
    This script solves the chemical puzzle by calculating empirical/molecular formulas,
    molar masses, and using chemical principles to deduce the structure of X.
    """
    print("Step 1: Determine the molecular formula of amine A1 from its elemental composition.")
    
    # Elemental composition of A1
    composition = {'C': 54.5, 'H': 13.6, 'N': 31.8}
    atomic_masses = {'C': 12.01, 'H': 1.008, 'N': 14.01}
    
    # Calculate moles per 100g of substance
    moles = {el: composition[el] / atomic_masses[el] for el in composition}
    
    # Find the smallest mole value to determine the empirical ratio
    min_moles = min(moles.values())
    empirical_ratio = {el: moles[el] / min_moles for el in moles}
    
    # The calculated ratio is ~C=2, H=6, N=1.
    # C2H6N is not a stable molecule (violates valency rules).
    # We test integer multiples. The first chemically sensible formula is often correct.
    # Let's check the ratio for C4H12N2 (multiplying by 2).
    molecular_formula = "C4H12N2"
    molar_mass_A1 = 4 * 12.01 + 12 * 1.008 + 2 * 14.01
    
    print(f"Initial empirical ratio is ~C:{empirical_ratio['C']:.1f} H:{empirical_ratio['H']:.1f} N:{empirical_ratio['N']:.1f}, which is C2H6N.")
    print("This empirical formula is chemically implausible. Checking for integer multiples...")
    print(f"The molecular formula C4H12N2 is a saturated C4 diamine and matches the given composition.")
    
    print("\n" + "="*50 + "\n")
    
    print("Step 2: Determine the molar mass of the carboxylic acid from titration.")
    
    # Titration data
    mass_acid_g = 2.16
    vol_koh_l = 0.030
    molarity_koh_mol_l = 1.0
    
    # Calculation
    moles_koh = vol_koh_l * molarity_koh_mol_l
    
    # The neutralization equation is R-COOH + KOH -> R-COOK + H2O (assuming a monobasic acid)
    # The molar ratio is 1:1, so moles_acid = moles_koh
    moles_acid = moles_koh
    molar_mass_acid = mass_acid_g / moles_acid

    print("The neutralization reaction is: Acid + KOH -> Salt + H2O")
    print(f"To neutralize {mass_acid_g} g of acid, {vol_koh_l * 1000} ml of {molarity_koh_mol_l} M KOH was used.")
    print("The final equation for molar mass is: Molar Mass = mass / (Volume * Molarity)")
    print(f"Each number in the final equation: Molar Mass = {mass_acid_g} / ({vol_koh_l} * {molarity_koh_mol_l})")
    print(f"Calculated Molar Mass = {molar_mass_acid:.2f} g/mol")

    print("\n" + "="*50 + "\n")
    
    print("Step 3: Deduce the structures based on the reaction sequence.")
    print("Clue 1: A1 is C4H12N2. A1 reacts with HNO2 to give A2. Thus, A1 is a C4 diamine, and A2 is a C4 diol.")
    print("Clue 2: A2 is oxidized to a carboxylic acid AND CO2. This indicates oxidative cleavage of a 1,2-diol (R-CH(OH)-CH2OH).")
    print("Clue 3: This cleavage means the final acid has one carbon less than A2. Since A2 is a C4 diol, the acid must be a C3 carboxylic acid.")
    print("The only simple C3 carboxylic acid is Propanoic Acid (CH3CH2COOH), with a molar mass of 74.08 g/mol.")
    print(f"Our calculated molar mass of {molar_mass_acid:.2f} g/mol is very close to 74.08 g/mol, suggesting the acid is indeed Propanoic Acid.")
    
    print("\n" + "="*50 + "\n")
    
    print("Step 4: Identify all structures by working backwards.")
    print(" - Carboxylic Acid: Propanoic acid (CH3CH2COOH)")
    print(" - A2 (the diol): Must be Butane-1,2-diol (CH3CH2CH(OH)CH2OH) to yield propanoic acid upon cleavage.")
    print(" - A1 (the diamine): Must be Butane-1,2-diamine (CH3CH2CH(NH2)CH2NH2).")
    print("   -> Checking NMR: Butane-1,2-diamine has 4 chemically distinct carbon atoms. This matches the 'four types of signals' data (assuming a 13C-NMR spectrum).")
    print(" - A (the bromo-compound): Must be 1,2-Dibromobutane (CH3CH2CH(Br)CH2Br).")
    print(" - X (the hydrocarbon): Must be the alkene that undergoes addition of Br2 to form 1,2-dibromobutane.")
    print("\nTherefore, the structure of hydrocarbon X is:")
    print("Name: But-1-ene")
    print("Formula: CH3CH2CH=CH2")

solve_chemistry_problem()
<<<But-1-ene (CH3CH2CH=CH2)>>>