import math

def reaction_calculator():
    """
    Calculates reagent quantities for the ethylation of 2-Methyl-1,4-naphthalenediol.
    """
    # --- Molar Masses (g/mol) ---
    MM_SM = 174.20       # 2-Methyl-1,4-naphthalenediol (C11H10O2)
    MM_NaH = 24.00       # Sodium Hydride (NaH)
    MM_EtBr = 108.97     # Ethyl bromide (C2H5Br)
    MM_Product = 230.31  # 1,4-Diethoxy-2-methylnaphthalene (C15H18O2)

    # --- Reaction Parameters from the user ---
    mass_sm = 10.0      # grams
    eq_NaH = 2.5        # equivalents of NaH
    eq_EtBr = 3.0       # equivalents of EtBr

    # --- Calculations ---
    # Moles of starting material (SM)
    moles_sm = mass_sm / MM_SM

    # Sodium Hydride (Base) calculations
    moles_NaH_to_add = moles_sm * eq_NaH
    mass_NaH_to_add = moles_NaH_to_add * MM_NaH
    # NaH is often a 60% dispersion in mineral oil
    mass_NaH_60_percent_to_add = mass_NaH_to_add / 0.60

    # Ethyl Bromide (Electrophile) calculations
    moles_EtBr_to_add = moles_sm * eq_EtBr
    mass_EtBr_to_add = moles_EtBr_to_add * MM_EtBr
    density_EtBr = 1.46  # g/mL
    volume_EtBr_to_add = mass_EtBr_to_add / density_EtBr

    # Theoretical product yield calculations
    # The reaction stoichiometry is 1 mol of SM to 1 mol of Product
    moles_product_theoretical = moles_sm
    mass_product_theoretical = moles_product_theoretical * MM_Product
    
    # --- Output Results ---
    print("--- Reaction Stoichiometry Plan ---")
    print(f"Starting Material: 2-Methyl-1,4-naphthalenediol")
    print(f"  - Mass: {mass_sm:.2f} g")
    print(f"  - Moles: {moles_sm:.4f} mol\n")

    print(f"Base: Sodium Hydride (NaH)")
    print(f"  - Equivalents: {eq_NaH}")
    print(f"  - Moles required: {moles_NaH_to_add:.4f} mol")
    print(f"  - Mass (100% NaH): {mass_NaH_to_add:.2f} g")
    print(f"  - Mass (60% dispersion in oil): {mass_NaH_60_percent_to_add:.2f} g\n")
    
    print(f"Electrophile: Ethyl Bromide (EtBr)")
    print(f"  - Equivalents: {eq_EtBr}")
    print(f"  - Moles required: {moles_EtBr_to_add:.4f} mol")
    print(f"  - Mass required: {mass_EtBr_to_add:.2f} g")
    print(f"  - Volume required: {volume_EtBr_to_add:.2f} mL\n")

    print(f"Product: 1,4-Diethoxy-2-methylnaphthalene")
    print(f"  - Theoretical Moles: {moles_product_theoretical:.4f} mol")
    print(f"  - Theoretical Mass (Yield): {mass_product_theoretical:.2f} g\n")
    
    # The balanced equation requires 2 equivalents of NaH and 2 of EtBr per SM.
    # C11H10O2 + 2 NaH + 2 EtBr -> C15H18O2 + 2 NaBr + 2 H2
    print("--- Final Equation with Reacting Masses ---")
    # We calculate the theoretical mass of EtBr that reacts (2 eq)
    mass_EtBr_reacted = moles_sm * 2 * MM_EtBr
    print(f"Equation: C11H10O2 + 2 EtBr -> C15H18O2 + 2 NaBr")
    print(f"Masses:   {mass_sm:.2f} g + {mass_EtBr_reacted:.2f} g -> {mass_product_theoretical:.2f} g")

if __name__ == '__main__':
    reaction_calculator()