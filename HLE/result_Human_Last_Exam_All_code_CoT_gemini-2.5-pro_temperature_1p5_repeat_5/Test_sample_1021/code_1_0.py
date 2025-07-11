def check_reaction_stoichiometry():
    """
    Calculates the molar quantities for the described SN2 reaction
    to verify if the stoichiometry was planned correctly.
    """
    # --- Given values ---
    mass_sm = 10.0  # grams of starting material
    eq_base_used = 2.5  # equivalents of NaH used
    eq_alkyl_used = 3.0  # equivalents of Ethyl Bromide used

    # --- Molecular Weights (g/mol) ---
    # Starting Material: 2-Methyl-1,4-naphthalenediol (C11H10O2)
    mw_sm = (11 * 12.011) + (10 * 1.008) + (2 * 15.999) # 174.20 g/mol

    # --- Calculations ---
    # Moles of Starting Material (SM)
    moles_sm = mass_sm / mw_sm

    # The reaction involves two hydroxyl groups, so it is a double SN2 reaction.
    # Stoichiometrically, 2 eq of base and 2 eq of alkylating agent are required.
    eq_base_needed = 2.0
    eq_alkyl_needed = 2.0

    # Moles of reagents used
    moles_base_used = moles_sm * eq_base_used
    moles_alkyl_used = moles_sm * eq_alkyl_used

    # --- Output Results ---
    print("--- Reaction Stoichiometry Analysis ---")
    print(f"Starting Material (SM): 2-Methyl-1,4-naphthalenediol")
    print(f"Molecular Weight of SM: {mw_sm:.2f} g/mol")
    print(f"Mass of SM used: {mass_sm} g")
    print(f"Moles of SM: {moles_sm:.4f} mol\n")

    print(f"Stoichiometric Requirement per mole of SM:")
    print(f"  - Base needed: {eq_base_needed} eq")
    print(f"  - Ethyl Bromide needed: {eq_alkyl_needed} eq\n")

    print(f"Actual Reagents Used:")
    print(f"  - Base (NaH) used: {eq_base_used} eq, which is {moles_base_used:.4f} mol")
    print(f"  - Ethyl Bromide used: {eq_alkyl_used} eq, which is {moles_alkyl_used:.4f} mol\n")

    print("Conclusion: The amounts of base and ethyl bromide used were in excess of the stoichiometric requirement.")
    print("This suggests the failure was not due to incorrect stoichiometry, pointing towards an issue like oxidation.")

check_reaction_stoichiometry()