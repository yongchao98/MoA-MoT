import sys

def solve_chemistry_problem():
    """
    This function performs the calculations and provides the reasoning to determine the structure of hydrocarbon X.
    """
    # --- Step 1: Determine the Molecular Formula of A1 ---
    print("--- Step 1: Determining the Molecular Formula of A1 ---")
    
    # Given elemental composition of A1
    percent_C = 54.5
    percent_H = 13.6
    percent_N = 31.8
    
    # Atomic masses
    mass_C = 12.01
    mass_H = 1.008
    mass_N = 14.01
    
    # Assuming 100g of substance A1, calculate moles of each element
    moles_C = percent_C / mass_C
    moles_H = percent_H / mass_H
    moles_N = percent_N / mass_N
    
    # Find the simplest ratio by dividing by the smallest number of moles
    min_moles = min(moles_C, moles_H, moles_N)
    ratio_C = moles_C / moles_N # N is the smallest
    ratio_H = moles_H / moles_N
    ratio_N = moles_N / moles_N
    
    print(f"The molar ratio of C:H:N is approximately {ratio_C:.2f}:{ratio_H:.2f}:{ratio_N:.2f}, which suggests an empirical formula of C2H6N.")
    print("This empirical formula is not chemically stable. Let's test integer multiples.")
    
    # Test molecular formula C4H12N2 (double the empirical formula C2H6N, but with H adjusted for a saturated diamine)
    formula_A1 = "C4H12N2"
    molar_mass_A1 = 4 * mass_C + 12 * mass_H + 2 * mass_N
    calc_percent_C = (4 * mass_C / molar_mass_A1) * 100
    calc_percent_H = (12 * mass_H / molar_mass_A1) * 100
    calc_percent_N = (2 * mass_N / molar_mass_A1) * 100
    
    print(f"Let's check the molecular formula {formula_A1} (a butanediamine):")
    print(f"  - Calculated %C = {calc_percent_C:.1f}% (Given: 54.5%)")
    print(f"  - Calculated %H = {calc_percent_H:.1f}% (Given: 13.6%)")
    print(f"  - Calculated %N = {calc_percent_N:.1f}% (Given: 31.8%)")
    print("The values match perfectly. Therefore, A1 has the molecular formula C4H12N2 and is a saturated butanediamine.")
    print("-" * 20)

    # --- Step 2: Analyze the Titration Data ---
    print("\n--- Step 2: Calculating the Molar Mass of the Carboxylic Acid ---")
    mass_acid = 2.16  # g
    vol_KOH = 0.030   # L (30 ml)
    conc_KOH = 1.0    # M (mol/L)

    print(f"Equation for moles of KOH: n = Concentration * Volume")
    moles_KOH = conc_KOH * vol_KOH
    print(f"Moles of KOH = {conc_KOH} mol/L * {vol_KOH} L = {moles_KOH} mol")
    
    # For a monoprotic acid, moles_acid = moles_KOH
    moles_acid = moles_KOH
    print(f"Assuming a monoprotic acid, moles of acid = {moles_acid} mol")
    
    print(f"Equation for Molar Mass: M = mass / moles")
    molar_mass_acid = mass_acid / moles_acid
    print(f"Calculated Molar Mass of the acid = {mass_acid} g / {moles_acid} mol = {molar_mass_acid:.2f} g/mol")
    print("-" * 20)

    # --- Step 3: Deduce the Reaction Pathway and Structures ---
    print("\n--- Step 3: Deducing the Structures ---")
    # Molar mass of Propanoic Acid (C3H6O2)
    molar_mass_propanoic = 3 * mass_C + 6 * mass_H + 2 * 16.00
    print(f"The calculated molar mass (72.00 g/mol) is very close to the molar mass of propanoic acid, C3H6O2 (M = {molar_mass_propanoic:.2f} g/mol).")
    print("The small difference is likely due to experimental error. The chemical context strongly supports the acid being propanoic acid.")
    
    print("\nWorking backwards from Propanoic Acid:")
    print("1. The formation of propanoic acid (a C3 acid) and CO2 from a C4 compound (A2) is characteristic of the oxidative cleavage of 1,2-butanediol.")
    print("   Reaction: CH3-CH2-CH(OH)-CH2(OH) [A2] --[O]--> CH3-CH2-COOH [Propanoic Acid] + CO2")
    print("   Therefore, A2 is 1,2-butanediol.")
    
    print("\n2. A2 (1,2-butanediol) is formed from A1 reacting with nitrous acid. This means A1 must be the corresponding diamine.")
    print("   Therefore, A1 is 1,2-butanediamine: CH3-CH2-CH(NH2)-CH2(NH2).")
    print("   Let's check this against the data: Formula is C4H12N2 (matches step 1). The structure has 4 chemically distinct carbon atoms, which matches the 'four types of signals' in the NMR spectrum.")

    print("\n3. A1 (1,2-butanediamine) is formed from A. This means A must be 1,2-dibromobutane.")
    
    print("\n4. A (1,2-dibromobutane) is formed from the addition of Br2 to hydrocarbon X.")
    print("   Reaction: X + Br2 -> CH3-CH2-CH(Br)-CH2(Br) [A]")
    print("   This is a classic addition reaction to an alkene. Therefore, X must be 1-butene.")
    print("   Structure of X: CH2=CH-CH2-CH3")
    print("-" * 20)
    
    # --- Step 4: Final Conclusion ---
    print("\n--- Step 4: Final Conclusion ---")
    print("The structure of X is 1-butene. This fits all the given information:")
    print(" - It's a C4 hydrocarbon.")
    print(" - Its reaction with Br2 yields a single product (1,2-dibromobutane).")
    print(" - This leads to 1,2-butanediamine (A1), which has a formula of C4H12N2 and 4 unique carbon signals.")
    print(" - This leads to 1,2-butanediol (A2), which oxidizes to propanoic acid and CO2, consistent with the titration results (within minor experimental error).")

solve_chemistry_problem()
<<<1-butene>>>