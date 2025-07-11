def calculate_reagent_volume():
    """
    Calculates the precise volume of n-BuLi needed for a reaction to avoid side products.

    The key to solving the problem of two boron signals is to use a precise stoichiometric
    amount (1.00 equivalents) of n-BuLi. This prevents the formation of excess aryllithium
    which can react with the desired product to form a diarylborinic acid byproduct.
    This function calculates the exact volume of n-BuLi solution required after its
    concentration has been accurately determined by titration.
    """

    # --- Reaction Parameters ---
    # Mass of the starting material (2-bromo-4-chloro-1-iodobenzene) in grams
    mass_sm = 5.0
    # Molecular weight of the starting material in g/mol
    mw_sm = 318.35
    # Precise concentration of n-BuLi solution from titration in mol/L (M)
    conc_nbuli = 1.55
    # Desired equivalents of n-BuLi to avoid side reactions
    equivalents_nbuli = 1.00

    # --- Calculation Steps ---
    # 1. Calculate moles of starting material (SM)
    moles_sm = mass_sm / mw_sm

    # 2. Calculate moles of n-BuLi needed for the desired stoichiometry
    moles_nbuli = moles_sm * equivalents_nbuli

    # 3. Calculate the required volume of n-BuLi solution in Liters
    volume_nbuli_L = moles_nbuli / conc_nbuli

    # 4. Convert volume to milliliters (mL) for practical use in the lab
    volume_nbuli_mL = volume_nbuli_L * 1000

    # --- Output the results with the final equation numbers ---
    print("To solve the problem, use a precise amount of n-BuLi (1.00 eq).")
    print("Here is the calculation to find the required volume:\n")
    
    print(f"1. Moles of Starting Material (SM) = Mass SM / MW SM")
    print(f"   Moles SM = {mass_sm:.2f} g / {mw_sm:.2f} g/mol")
    print(f"   Moles SM = {moles_sm:.5f} mol\n")

    print(f"2. Moles of n-BuLi required = Moles SM * Equivalents")
    print(f"   Moles n-BuLi = {moles_sm:.5f} mol * {equivalents_nbuli:.2f}")
    print(f"   Moles n-BuLi = {moles_nbuli:.5f} mol\n")

    print(f"3. Volume of n-BuLi = Moles n-BuLi / Concentration n-BuLi")
    print(f"   Volume (L) = {moles_nbuli:.5f} mol / {conc_nbuli:.2f} M")
    print(f"   Volume (L) = {volume_nbuli_L:.5f} L\n")

    print("--------------------------------------------------")
    print(f"Final Precise Volume to add: {volume_nbuli_mL:.2f} mL")
    print("--------------------------------------------------")
    print("\nUsing this precisely calculated volume minimizes the formation of the diaryl byproduct,")
    print("which is the cause of the second boron signal in the NMR.")

# Execute the calculation
calculate_reagent_volume()