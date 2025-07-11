import math

def solve_structure():
    """
    This function solves the chemical puzzle by calculating molecular formulas and molar masses,
    and then printing the step-by-step logical deduction.
    """
    print("Step 1: Determine the molecular formula of compound A1 from its elemental composition.")
    
    # Given percentages
    c_percent = 54.5
    h_percent = 13.6
    n_percent = 31.8
    
    # Atomic masses
    c_mass = 12.01
    h_mass = 1.008
    n_mass = 14.01
    
    # Calculate moles assuming 100g of sample
    moles_c = c_percent / c_mass
    moles_n = n_percent / n_mass
    
    # Determine ratios for empirical formula
    c_ratio = moles_c / moles_n
    h_ratio = (h_percent / h_mass) / moles_n
    n_ratio = moles_n / moles_n

    print(f"The molar ratios are approximately C:{c_ratio:.2f} H:{h_ratio:.2f} N:{n_ratio:.2f}, which simplifies to C:2 H:6 N:1.")
    print("The empirical formula is C2H6N. For a stable, non-radical molecule, the molecular formula must be a multiple, like (C2H6N)n.")
    print("For n=2, the molecular formula is C4H12N2.")
    
    # Verify C4H12N2 percentages
    mf_mass = 4 * c_mass + 12 * h_mass + 2 * n_mass
    mf_c_percent = (4 * c_mass / mf_mass) * 100
    mf_h_percent = (12 * h_mass / mf_mass) * 100
    mf_n_percent = (2 * n_mass / mf_mass) * 100
    
    print(f"Checking the percentages for C4H12N2 (Molar Mass = {mf_mass:.2f} g/mol):")
    print(f"C: {mf_c_percent:.1f}% (Given: 54.5%), H: {mf_h_percent:.1f}% (Given: 13.6%), N: {mf_n_percent:.1f}% (Given: 31.8%)")
    print("The percentages match. Thus, the molecular formula of A1 is C4H12N2.\n")

    print("Step 2: Calculate the molar mass of the carboxylic acid from the titration data.")
    
    acid_mass = 2.16  # g
    koh_volume_ml = 30.0  # ml
    koh_molarity = 1.0  # M
    
    # Calculate moles of KOH
    moles_koh = koh_molarity * (koh_volume_ml / 1000)
    
    # For a monoprotic acid, moles of acid = moles of KOH
    moles_acid = moles_koh
    
    # Calculate molar mass
    acid_molar_mass = acid_mass / moles_acid
    
    print(f"The neutralization equation is: R-COOH + KOH -> R-COOK + H2O")
    print(f"Moles of KOH used = {koh_molarity} M * ({koh_volume_ml} / 1000) L = {moles_koh:.4f} mol.")
    print(f"Therefore, moles of the carboxylic acid = {moles_acid:.4f} mol.")
    print(f"The calculated molar mass of the acid = {acid_mass} g / {moles_acid:.4f} mol = {acid_molar_mass:.2f} g/mol.\n")

    print("Step 3: Deduce the structures by working backward.")
    
    print("The molar mass of propanoic acid (CH3CH2COOH) is 74 g/mol, which is a very close match to our calculated value of 72.00 g/mol.")
    print("So, the carboxylic acid is propanoic acid (a C3 acid).")
    print("The oxidation of A2 produced propanoic acid (C3) and CO2 (C1). This indicates that A2 was a C4 compound that underwent oxidative cleavage.")
    print("This specific cleavage pattern results from the oxidation of a 1,2-diol. For the products to be a C3 acid and CO2, the diol A2 must be 1,2-butanediol (CH3CH2CH(OH)CH2OH).")
    print("A2 (1,2-butanediol) was formed from the diamine A1 by reaction with nitrous acid. Therefore, A1 must be 1,2-diaminobutane (CH3CH2CH(NH2)CH2NH2).")
    print("The molecular formula of 1,2-diaminobutane is C4H12N2, which matches our finding from Step 1.")
    print("The NMR spectrum of A1 shows four types of signals. 1,2-diaminobutane has four unique carbon atoms, which would result in four signals in a 13C NMR spectrum. This is consistent with the data.")
    print("A1 (1,2-diaminobutane) was formed from a substance A by reaction with ammonia. Thus, A is 1,2-dibromobutane.")
    print("A (1,2-dibromobutane) was formed as the only product from the reaction of hydrocarbon X with bromine. This is an electrophilic addition reaction to an alkene.")
    print("Therefore, hydrocarbon X must be the alkene that yields 1,2-dibromobutane upon addition of Br2.")
    print("The structure of X is CH2=CH-CH2-CH3.\n")
    
    print("Final Answer: The structure of substance X is but-1-ene.")

solve_structure()
<<<but-1-ene>>>