def calculate_reagent_quantities():
    """
    Calculates the molar and mass amounts of reagents for the described SN2 reaction.
    """
    # --- Given Information ---
    mass_sm = 10.0  # grams of starting material
    eq_nah = 2.5    # equivalents of Sodium Hydride
    eq_etbr = 3.0   # equivalents of Ethyl Bromide

    # --- Molecular Weights (g/mol) ---
    mw_sm = 174.20   # 2-Methyl-1,4-naphthalenediol (C11H10O2)
    mw_nah = 24.00   # Sodium Hydride (NaH)
    mw_etbr = 108.97 # Ethyl Bromide (C2H5Br)

    # --- Calculations ---
    # 1. Moles of Starting Material (SM)
    moles_sm = mass_sm / mw_sm

    # 2. Moles and Mass of Sodium Hydride (NaH)
    moles_nah = moles_sm * eq_nah
    mass_nah = moles_nah * mw_nah

    # 3. Moles and Mass of Ethyl Bromide (EtBr)
    moles_etbr = moles_sm * eq_etbr
    mass_etbr = moles_etbr * mw_etbr

    # --- Print Results ---
    print("--- Reaction Stoichiometry ---")
    print(f"Starting Material: 2-Methyl-1,4-naphthalenediol")
    print(f"Molecular Weight: {mw_sm} g/mol")
    print(f"Starting Mass: {mass_sm:.2f} g")
    print(f"Amount: {moles_sm:.4f} mol\n")

    print(f"Base: Sodium Hydride (NaH)")
    print(f"Molecular Weight: {mw_nah} g/mol")
    print(f"Equivalents used: {eq_nah}")
    print(f"Required Amount: {moles_nah:.4f} mol")
    print(f"Required Mass: {mass_nah:.2f} g\n")

    print(f"Electrophile: Ethyl Bromide (EtBr)")
    print(f"Molecular Weight: {mw_etbr} g/mol")
    print(f"Equivalents used: {eq_etbr}")
    print(f"Required Amount: {moles_etbr:.4f} mol")
    print(f"Required Mass: {mass_etbr:.2f} g")

calculate_reagent_quantities()