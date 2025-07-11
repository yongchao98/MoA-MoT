import sys

def solve_chemistry_problem():
    """
    This script solves the multi-step chemistry problem to identify hydrocarbon X.
    It follows a logical deduction process, starting from the final products and working backwards.
    """
    # Use simple integer atomic masses as is common in such problems
    # and to align with the likely calculation method for the given data.
    atomic_mass = {'C': 12, 'H': 1, 'N': 14, 'O': 16, 'K': 39, 'Br': 80}

    print("Step 1: Determine the molecular formula of substance A1.")
    print("Given composition: C=54.5%, H=13.6%, N=31.8%")
    print("Assuming a 100g sample, we calculate the moles of each element:")
    
    comp_C, comp_H, comp_N = 54.5, 13.6, 31.8
    mass_C, mass_H, mass_N = atomic_mass['C'], atomic_mass['H'], atomic_mass['N']

    moles_C = comp_C / mass_C
    moles_H = comp_H / mass_H
    moles_N = comp_N / mass_N
    print(f"Moles C = {comp_C} / {mass_C} = {moles_C:.3f}")
    print(f"Moles H = {comp_H} / {mass_H} = {moles_H:.3f}")
    print(f"Moles N = {comp_N} / {mass_N} = {moles_N:.3f}")

    min_moles = min(moles_C, moles_H, moles_N)
    print(f"\nDividing by the smallest value ({min_moles:.3f}) to find the simplest ratio:")
    ratio_C = moles_C / min_moles
    ratio_H = moles_H / min_moles
    ratio_N = moles_N / min_moles
    print(f"C: {ratio_C:.1f}, H: {ratio_H:.1f}, N: {ratio_N:.1f}")
    print("The empirical formula is C2H6N.")

    print("\nA1 is a stable organic molecule formed from a reaction with ammonia, suggesting it is an amine.")
    print("An even number of nitrogen atoms is typical for stable amines formed this way. The simplest molecular formula is (C2H6N)x2 = C4H12N2 (a diamine).")
    mass_C4H12N2 = 4*mass_C + 12*mass_H + 2*mass_N
    check_C = (4*mass_C/mass_C4H12N2)*100
    check_H = (12*mass_H/mass_C4H12N2)*100
    check_N = (2*mass_N/mass_C4H12N2)*100
    print(f"Checking the percentages for C4H12N2 (Molar Mass = {mass_C4H12N2}): C={check_C:.1f}%, H={check_H:.1f}%, N={check_N:.1f}%. This matches the given data.")
    print("Therefore, the molecular formula of A1 is C4H12N2.\n")


    print("Step 2: Calculate the molar mass of the carboxylic acid.")
    mass_acid = 2.16  # g
    vol_KOH_L = 0.030 # L (30 ml)
    conc_KOH = 1.0    # M
    
    print("The reaction is a 1:1 neutralization: Acid + KOH -> Salt + H2O.")
    moles_KOH = vol_KOH_L * conc_KOH
    print(f"Moles of KOH used = {vol_KOH_L:.3f} L * {conc_KOH:.1f} M = {moles_KOH:.3f} mol.")
    
    # Assuming it is a monoprotic acid
    moles_acid = moles_KOH
    print(f"This means moles of carboxylic acid = {moles_acid:.3f} mol.")
    
    molar_mass_acid = mass_acid / moles_acid
    print(f"The calculated molar mass of the acid is: {mass_acid} g / {moles_acid:.2f} mol = {molar_mass_acid:.0f} g/mol.\n")


    print("Step 3: Deduce the structures of the carboxylic acid and substance A2.")
    propanoic_acid_mass = 3 * mass_C + 6 * mass_H + 2 * mass_O
    print(f"The calculated molar mass (72 g/mol) is very close to that of propanoic acid, C3H6O2 (Molar Mass = 3*{mass_C} + 6*{mass_H} + 2*{atomic_mass['O']} = {propanoic_acid_mass} g/mol).")
    print("The formation of a C3 acid (propanoic acid) and CO2 from a C4 precursor (A2) strongly suggests the oxidative cleavage of a vicinal diol.")
    print("Specifically, the strong oxidation of 1,2-butanediol (CH3CH2CH(OH)CH2(OH)) yields propanoic acid (CH3CH2COOH) and CO2 (from the -CH2OH carbon).")
    print("Therefore, the carboxylic acid is Propanoic Acid and substance A2 is 1,2-Butanediol.\n")
    

    print("Step 4: Deduce the structure of A1.")
    print("A2 (1,2-butanediol) is formed from A1 by reaction with nitrous acid, which converts primary amino groups (-NH2) to hydroxyl groups (-OH).")
    print("Thus, A1 must be 1,2-butanediamine: CH3CH2CH(NH2)CH2(NH2).")
    print("This structure corresponds to the molecular formula C4H12N2 and has four distinct sets of C-H protons (from CH3, CH2, CH, and CH2NH2), which matches the 'four types of signals' in the NMR spectrum.\n")
    

    print("Step 5: Deduce the structure of A.")
    print("A1 (1,2-butanediamine) is formed from substance A by reaction with excess ammonia, where bromine atoms are replaced by amino groups.")
    print("Therefore, A must be the corresponding dibromoalkane: 1,2-dibromobutane, CH3CH2CH(Br)CH2(Br).\n")


    print("Step 6: Deduce the structure of X.")
    print("A (1,2-dibromobutane) is the sole product of the reaction between hydrocarbon X and bromine (Br2).")
    print("This is a classic electrophilic addition of bromine across a double bond.")
    print("Therefore, the original hydrocarbon X must be But-1-ene: CH3CH2CH=CH2.\n")
    
    print("Final answer for the structure of X:")


solve_chemistry_problem()
sys.stdout.flush()
print('<<<But-1-ene>>>')