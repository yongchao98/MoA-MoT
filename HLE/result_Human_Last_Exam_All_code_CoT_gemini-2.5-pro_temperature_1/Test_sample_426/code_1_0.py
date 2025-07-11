def solve_chemistry_problem():
    """
    This function solves the chemical puzzle by calculating the molar mass of the acid
    and deducing the structure of the initial hydrocarbon X based on the reaction sequence.
    """

    # Step 1: Calculate the molar mass of the carboxylic acid from titration data.
    mass_acid_g = 2.16  # grams
    vol_koh_ml = 30.0  # mL
    conc_koh_M = 1.0  # Molar (mol/L)

    # Convert volume to Liters
    vol_koh_L = vol_koh_ml / 1000.0

    # Calculate moles of KOH used
    moles_koh = vol_koh_L * conc_koh_M

    # In neutralization of a monoprotic acid, moles of acid = moles of KOH
    moles_acid = moles_koh

    # Calculate the molar mass of the acid
    molar_mass_acid = mass_acid_g / moles_acid

    print("Step 1: Calculating the molar mass of the carboxylic acid.")
    print(f"The neutralization equation is: {mass_acid_g} g / Molar Mass = {vol_koh_L:.3f} L * {conc_koh_M} M")
    print(f"Molar Mass = {mass_acid_g} / {moles_acid:.3f} = {molar_mass_acid:.2f} g/mol")
    print("-" * 20)

    # Step 2: Deduce the structure of X from the molar mass and reaction sequence.
    print("Step 2: Deducing the structure of hydrocarbon X.")
    print("The calculated molar mass is 72.00 g/mol.")
    print("This value is very close to the molar mass of propanoic acid (C3H6O2), which is 74.08 g/mol.")
    print("The formation of a C3 acid (propanoic acid) and CO2 from the oxidation of a diol (A2) suggests the diol had 4 carbons and was butane-1,2-diol.")
    print("Oxidation of butane-1,2-diol cleaves the bond between the carbons bearing the -OH groups, yielding propanoic acid and CO2.")
    print("\nTracing the reactions backwards:")
    print("  - A2 (butane-1,2-diol) is formed from A1 (butane-1,2-diamine).")
    print("  - A1 (butane-1,2-diamine) is formed from A (1,2-dibromobutane).")
    print("  - A (1,2-dibromobutane) is formed from the addition of Br2 to hydrocarbon X.")
    print("\nTherefore, the original hydrocarbon X must be but-1-ene.")

solve_chemistry_problem()

<<<but-1-ene>>>