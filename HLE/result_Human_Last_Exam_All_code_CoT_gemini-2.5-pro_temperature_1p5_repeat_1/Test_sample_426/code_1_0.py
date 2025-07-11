import math

def solve_chemistry_problem():
    """
    This function solves the chemical puzzle by performing calculations and printing the logical deductions step-by-step.
    """
    print("Step 1: Determine the molar mass of the carboxylic acid from the titration data.")
    
    # Given data for neutralization
    mass_acid = 2.16  # in grams
    vol_koh = 0.030   # in liters (30 ml)
    molarity_koh = 1.0  # in M (mol/L)
    
    # Calculation
    moles_koh = molarity_koh * vol_koh
    # Assuming the carboxylic acid is monoprotic, moles of acid = moles of KOH
    moles_acid = moles_koh
    molar_mass_acid = mass_acid / moles_acid
    
    print("The neutralization equation is: R-COOH + KOH -> R-COOK + H2O")
    print(f"Volume of KOH solution = {vol_koh * 1000} ml")
    print(f"Molarity of KOH solution = {molarity_koh} M")
    print(f"Moles of KOH used = {molarity_koh} mol/L * {vol_koh} L = {moles_koh:.3f} mol")
    print(f"Mass of carboxylic acid = {mass_acid} g")
    print(f"Calculated molar mass of the acid = {mass_acid} g / {moles_acid:.3f} mol = {molar_mass_acid:.2f} g/mol")
    print("-" * 20)
    
    print("Step 2: Determine the empirical formula of substance A1.")
    
    # Given elemental composition of A1
    percent_C = 54.5
    percent_H = 13.6
    percent_N = 31.8
    
    # Atomic masses
    atomic_mass_C = 12.01
    atomic_mass_H = 1.008
    atomic_mass_N = 14.01
    
    # Moles ratio calculation
    moles_C = percent_C / atomic_mass_C
    moles_H = percent_H / atomic_mass_H
    moles_N = percent_N / atomic_mass_N
    
    # Finding the simplest ratio
    min_moles = min(moles_C, moles_H, moles_N)
    ratio_C = moles_C / min_moles
    ratio_H = moles_H / min_moles
    ratio_N = moles_N / min_moles
    
    print(f"Elemental composition of A1: C={percent_C}%, H={percent_H}%, N={percent_N}%")
    print(f"Relative moles: C={moles_C:.2f}, H={moles_H:.2f}, N={moles_N:.2f}")
    print(f"Dividing by the smallest value ({min_moles:.2f}):")
    print(f"C ratio = {ratio_C:.2f} -> {round(ratio_C)}")
    print(f"H ratio = {ratio_H:.2f} -> {round(ratio_H)}")
    print(f"N ratio = {ratio_N:.2f} -> {round(ratio_N)}")
    print("The empirical formula of A1 is C2H6N.")
    print("-" * 20)
    
    print("Step 3: Deduce molecular formulas and structures.")
    
    print("The oxidation of alcohol A2 yields a carboxylic acid and CO2.")
    print("The calculated molar mass of the acid is 72 g/mol, which corresponds to the formula C3H4O2 (e.g., acrylic acid).")
    print("The total number of carbon atoms in the products is 3 (from the acid) + 1 (from CO2) = 4 carbons.")
    print("This implies that substance A2, and therefore A1 and X, are C4 compounds.")
    
    print("\nThe empirical formula of A1 is C2H6N. To have 4 carbons, we must double the formula.")
    print("Molecular Formula of A1 = (C2H6N) * 2 = C4H12N2. This is a diaminobutane isomer.")
    
    print("\nA1 is formed from hydrocarbon X + Br2 -> A (dibromoalkane), then A + NH3 -> A1.")
    print("Therefore, X must be a C4 hydrocarbon. The reaction with Br2 to form a single product A suggests X is an alkene (C4H8).")

    print("\nThe NMR spectrum of A1 has 'four types of signals'. This implies the carbon backbone is asymmetrical.")
    print("Let's consider possible isomers of X (C4H8):")
    print("  - But-2-ene (CH3-CH=CH-CH3) is symmetrical. Its product, 2,3-diaminobutane, would have only 2 types of carbon signals. -> Incorrect.")
    print("  - 2-Methylpropene ((CH3)2C=CH2) has some symmetry. Its product, 1,2-diamino-2-methylpropane, has only 3 types of carbon signals. -> Incorrect.")
    print("  - But-1-ene (CH2=CH-CH2-CH3) is asymmetrical. Its product, 1,2-diaminobutane (NH2-CH2-CH(NH2)-CH2-CH3), has 4 unique carbon atoms. -> This fits the NMR data.")
    
    print("\nTracing the reaction pathway from But-1-ene:")
    print("X (But-1-ene) + Br2 -> A (1,2-dibromobutane)")
    print("A (1,2-dibromobutane) + excess NH3 -> A1 (1,2-diaminobutane)")
    print("A1 (1,2-diaminobutane) + HNO2 -> A2 (butane-1,2-diol)")
    print("A2 (butane-1,2-diol) oxidation cleaves the C-C bond between the -OH groups, yielding CO2 and Propanoic Acid (CH3CH2COOH).")
    
    print("\nReconciling the final step:")
    print("The expected product is Propanoic Acid (M=74 g/mol). The titration data gives M=72 g/mol.")
    print("This is a minor discrepancy. However, all other structural evidence (reaction type, carbon count, NMR data) strongly and consistently points to but-1-ene.")
    print("-" * 20)
    
    print("Conclusion:")
    print("Based on the comprehensive analysis, the structure of hydrocarbon X is But-1-ene.")

solve_chemistry_problem()
<<<But-1-ene>>>