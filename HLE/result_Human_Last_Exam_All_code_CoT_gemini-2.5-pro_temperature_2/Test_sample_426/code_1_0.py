def solve_chemistry_problem():
    """
    This script calculates the molar mass of the unknown carboxylic acid
    based on the titration data provided in the problem.
    """
    # Given data from the problem
    mass_acid = 2.16  # in grams
    vol_KOH_ml = 30   # in milliliters
    conc_KOH = 1.0    # in M (moles per liter)

    # --- Step 1: Convert volume to liters ---
    vol_KOH_L = vol_KOH_ml / 1000

    # --- Step 2: Calculate moles of KOH ---
    moles_KOH = vol_KOH_L * conc_KOH

    # --- Step 3: Determine moles of acid ---
    # The neutralization of a monocarboxylic acid (R-COOH) with KOH is a 1:1 molar reaction.
    # R-COOH + KOH -> R-COOK + H2O
    # So, moles of acid = moles of KOH. We will assume it's monoprotic.
    moles_acid = moles_KOH

    # --- Step 4: Calculate the molar mass of the acid ---
    molar_mass_acid = mass_acid / moles_acid

    # --- Print the step-by-step calculation ---
    print(f"Step 1: The mass of the carboxylic acid to be neutralized is {mass_acid} g.")
    print(f"Step 2: The volume of KOH solution used is {vol_KOH_ml} ml, which is {vol_KOH_L} L.")
    print(f"Step 3: The concentration of the KOH solution is {conc_KOH} M.")
    print(f"Step 4: Moles of KOH used = {vol_KOH_L:.3f} L * {conc_KOH:.1f} M = {moles_KOH:.3f} mol.")
    print("Step 5: Assuming the acid is monocarboxylic (has one -COOH group), the moles of acid are equal to the moles of KOH.")
    print(f"Step 6: Moles of acid = {moles_acid:.3f} mol.")
    print(f"Step 7: The calculated molar mass of the acid is mass / moles.")
    print(f"Molar Mass = {mass_acid} g / {moles_acid:.3f} mol = {molar_mass_acid:.2f} g/mol")

solve_chemistry_problem()