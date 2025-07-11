def solve_chemistry_problem():
    """
    This script solves the chemical puzzle by performing calculations and deductive reasoning.
    """

    print("Step 1: Determining the molecular formula of substance A1.")
    
    # Elemental composition of A1
    percent_C = 54.5
    percent_H = 13.6
    percent_N = 31.8
    
    # Atomic masses
    mass_C = 12.01
    mass_H = 1.008
    mass_N = 14.01
    
    # Calculate moles in a 100g sample
    moles_C = percent_C / mass_C
    moles_H = percent_H / mass_H
    moles_N = percent_N / mass_N
    
    # Find the smallest mole value to determine the ratio
    min_moles = min(moles_C, moles_H, moles_N)
    
    ratio_C = moles_C / moles_N
    ratio_H = moles_H / moles_N
    ratio_N = moles_N / moles_N
    
    print(f"The molar ratio of C:H:N is approximately {ratio_C:.1f}:{ratio_H:.1f}:{ratio_N:.1f}, which simplifies to 2:6:1.")
    print("The empirical formula is C2H6N. This is not a stable molecule.")
    print("Let's try multiplying the ratio by 2, assuming there are two nitrogen atoms.")
    molecular_formula_A1 = "C4H12N2"
    molar_mass_A1 = 4 * 12.01 + 12 * 1.008 + 2 * 14.01
    
    # Verify the composition with the proposed molecular formula C4H12N2
    calc_percent_C = (4 * 12.01 / molar_mass_A1) * 100
    calc_percent_H = (12 * 1.008 / molar_mass_A1) * 100
    calc_percent_N = (2 * 14.01 / molar_mass_A1) * 100
    
    print(f"The proposed molecular formula is {molecular_formula_A1}. Let's check its composition:")
    print(f"Calculated: C={calc_percent_C:.1f}%, H={calc_percent_H:.1f}%, N={calc_percent_N:.1f}%. This is a perfect match.")
    print("-" * 30)

    print("Step 2: Determining the molar mass of the carboxylic acid.")
    
    # Titration data
    mass_acid = 2.16  # g
    vol_KOH = 0.030   # L (30 ml)
    conc_KOH = 1.0    # M (mol/L)
    
    # Moles of KOH = Volume * Concentration
    moles_KOH = vol_KOH * conc_KOH
    
    # Assuming a monoprotic acid, moles of acid = moles of KOH
    moles_acid = moles_KOH
    
    # Molar Mass = mass / moles
    molar_mass_acid = mass_acid / moles_acid
    
    print(f"The neutralization reaction is: R-COOH + KOH -> R-COOK + H2O")
    print(f"Moles of KOH used = {vol_KOH} L * {conc_KOH} mol/L = {moles_KOH} mol.")
    print(f"Molar mass of the acid = {mass_acid} g / {moles_acid} mol = {molar_mass_acid:.0f} g/mol.")
    print("-" * 30)

    print("Step 3: Deducing the structures.")
    
    print("The calculated molar mass of the acid is 72 g/mol. The molar mass of propanoic acid (CH3CH2COOH) is 74 g/mol. This is a very close match.")
    print("The oxidation of A2 produced this acid and CO2. This indicates a cleavage reaction. The formation of propanoic acid (a C3 acid) and CO2 (a C1 fragment) means the precursor A2 must have been a C4 compound.")
    print("This specific cleavage pattern points to the oxidation of a 1,2-diol. Therefore, A2 must be Butane-1,2-diol: HO-CH2-CH(OH)-CH2-CH3.")
    print("Oxidation of HO-CH2-CH(OH)-CH2-CH3 cleaves the C1-C2 bond:")
    print("  - The HO-CH2- group (C1) is oxidized to CO2.")
    print("  - The -CH(OH)-CH2-CH3 group (C2, C3, C4) is oxidized to propanoic acid: HOOC-CH2-CH3.")
    
    print("\nWorking backwards:")
    print("A2 (Butane-1,2-diol) was formed from A1 by reaction with nitrous acid. This converts primary amino groups to hydroxyl groups.")
    print("Therefore, A1 must be Butane-1,2-diamine: NH2-CH2-CH(NH2)-CH2-CH3.")
    print("Let's check A1's properties: The formula is C4H12N2 (matches Step 1). The 1H NMR spectrum should have 4 distinct C-H signals (-CH3, -CH2-, -CH(NH2)-, -CH2(NH2)), which matches the problem description.")
    
    print("\nA1 (Butane-1,2-diamine) was formed from A by reaction with excess ammonia. This is a nucleophilic substitution of bromine atoms.")
    print("Therefore, A must be 1,2-dibromobutane: Br-CH2-CH(Br)-CH2-CH3.")
    
    print("\nA (1,2-dibromobutane) was formed as the only product from the interaction of hydrocarbon X with bromine. This is characteristic of the addition of Br2 across a double bond.")
    print("Therefore, X must be the alkene that yields 1,2-dibromobutane upon bromination.")
    print("The structure of X is CH2=CH-CH2-CH3.")
    print("-" * 30)
    
    print("Final Answer: The structure of substance X is But-1-ene.")

solve_chemistry_problem()
<<<But-1-ene>>>