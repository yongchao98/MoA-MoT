def calculate_reagents():
    """
    Calculates the stoichiometry for the ethylation of 2-Methyl-1,4-naphthalenediol.
    This script helps in preparing for the reaction by calculating the required
    mass, moles, and volume of each component.
    """

    # --- Molecular Weights (g/mol) ---
    mw_sm = 174.20  # Starting Material: 2-Methyl-1,4-naphthalenediol (C11H10O2)
    mw_nah = 24.00   # Sodium Hydride (NaH)
    mw_etbr = 108.97 # Ethyl Bromide (EtBr)
    mw_product = 230.30 # Product: 2-Methyl-1,4-diethoxynaphthalene (C15H18O2)

    # --- Densities (g/mL) ---
    density_etbr = 1.46

    # --- Reaction Scale and Equivalents ---
    mass_sm = 10.0  # grams of starting material
    eq_nah = 2.5    # equivalents of NaH
    eq_etbr = 3.0   # equivalents of EtBr

    # --- Calculations ---
    # 1. Moles of Starting Material
    moles_sm = mass_sm / mw_sm

    # 2. Sodium Hydride (NaH)
    moles_nah = moles_sm * eq_nah
    mass_nah_pure = moles_nah * mw_nah
    # NaH is often supplied as a 60% dispersion in mineral oil
    mass_nah_60_percent = mass_nah_pure / 0.60

    # 3. Ethyl Bromide (EtBr)
    moles_etbr = moles_sm * eq_etbr
    mass_etbr = moles_etbr * mw_etbr
    volume_etbr = mass_etbr / density_etbr
    
    # 4. Theoretical Yield of Product
    # The reaction is 1:1 mole ratio between starting material and product
    moles_product_theoretical = moles_sm
    mass_product_theoretical = moles_product_theoretical * mw_product

    # --- Output Results ---
    print("--- Reaction Stoichiometry Calculation ---")
    print("\n[Starting Material: 2-Methyl-1,4-naphthalenediol]")
    print(f"Mass: {mass_sm:.2f} g")
    # Equation: moles = mass / molar mass
    print(f"Calculation for Moles: {mass_sm:.2f} g / {mw_sm:.2f} g/mol = {moles_sm:.4f} moles")

    print("\n[Base: Sodium Hydride (NaH)]")
    print(f"Equivalents: {eq_nah}")
    # Equation: moles_reagent = moles_sm * equivalents
    print(f"Calculation for Moles: {moles_sm:.4f} moles * {eq_nah:.1f} eq = {moles_nah:.4f} moles")
    print(f"Required Mass (Pure NaH): {mass_nah_pure:.2f} g")
    print(f"Required Mass (60% dispersion in oil): {mass_nah_60_percent:.2f} g")

    print("\n[Electrophile: Ethyl Bromide (EtBr)]")
    print(f"Equivalents: {eq_etbr}")
    # Equation: moles_reagent = moles_sm * equivalents
    print(f"Calculation for Moles: {moles_sm:.4f} moles * {eq_etbr:.1f} eq = {moles_etbr:.4f} moles")
    print(f"Required Mass: {mass_etbr:.2f} g")
    print(f"Required Volume: {volume_etbr:.2f} mL")
    
    print("\n[Product: 2-Methyl-1,4-diethoxynaphthalene]")
    print(f"Theoretical Yield (Mass): {mass_product_theoretical:.2f} g")

if __name__ == '__main__':
    calculate_reagents()