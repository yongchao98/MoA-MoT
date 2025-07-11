import math

def calculate_reaction_stoichiometry():
    """
    Calculates the required mass and moles for the ethylation of 2-Methyl-1,4-naphthalenediol.
    """
    # --- Constants (Molecular Weights in g/mol and Density) ---
    MW_STARTING_MATERIAL = (11 * 12.011) + (10 * 1.008) + (2 * 15.999)  # C11H10O2
    MW_NAH = 22.990 + 1.008  # NaH
    MW_ETBR = (2 * 12.011) + (5 * 1.008) + 79.904  # C2H5Br
    DENSITY_ETBR = 1.46  # g/mL
    MW_PRODUCT = (15 * 12.011) + (18 * 1.008) + (2 * 15.999) # C15H18O2

    # --- Inputs ---
    mass_starting_material_g = 10.0
    eq_nah = 2.5
    eq_etbr = 3.0
    
    # --- Calculations ---
    # 1. Moles of Starting Material (SM)
    moles_sm = mass_starting_material_g / MW_STARTING_MATERIAL
    
    # 2. Sodium Hydride (NaH) - Base
    moles_nah = moles_sm * eq_nah
    # Note: Commercial NaH is often a 60% dispersion in mineral oil.
    mass_nah_60_percent_dispersion_g = (moles_nah * MW_NAH) / 0.60
    
    # 3. Ethyl Bromide (EtBr) - Electrophile
    moles_etbr = moles_sm * eq_etbr
    mass_etbr_g = moles_etbr * MW_ETBR
    volume_etbr_ml = mass_etbr_g / DENSITY_ETBR
    
    # 4. Product - Theoretical Yield
    # The reaction is 1:1 between starting material and product
    moles_product_theoretical = moles_sm
    mass_product_theoretical_g = moles_product_theoretical * MW_PRODUCT
    
    # --- Output Results ---
    print("--- Reaction Stoichiometry Calculation ---")
    print(f"Starting Material (2-Methyl-1,4-naphthalenediol, MW={MW_STARTING_MATERIAL:.2f} g/mol):")
    print(f"  Mass: {mass_starting_material_g:.2f} g")
    print(f"  Moles: {moles_sm:.4f} mol")
    print("-" * 30)
    
    print(f"Sodium Hydride (NaH, {eq_nah} eq):")
    print(f"  Required Moles: {moles_nah:.4f} mol")
    print(f"  Required Mass (60% dispersion): {mass_nah_60_percent_dispersion_g:.2f} g")
    print("-" * 30)

    print(f"Ethyl Bromide (EtBr, {eq_etbr} eq, MW={MW_ETBR:.2f} g/mol):")
    print(f"  Required Moles: {moles_etbr:.4f} mol")
    print(f"  Required Mass: {mass_etbr_g:.2f} g")
    print(f"  Required Volume: {volume_etbr_ml:.2f} mL")
    print("-" * 30)

    print("--- Balanced Equation with Molar Equivalents ---")
    # The equation shows the molar equivalents relative to the starting material.
    # The reaction consumes 2 equivalents of base and 2 equivalents of electrophile for full conversion.
    # The student used 2.5 eq of NaH and 3.0 eq of EtBr.
    print(f"1.0 C11H10O2 + {eq_nah:.1f} NaH + {eq_etbr:.1f} EtBr -> 1.0 C15H18O2 + 2 NaBr + 0.5 NaH + 1.0 EtBr + H2")
    print("-" * 30)

    print("--- Theoretical Yield ---")
    print(f"Product (1,4-Diethoxy-2-methylnaphthalene, MW={MW_PRODUCT:.2f} g/mol):")
    print(f"  Theoretical Moles: {moles_product_theoretical:.4f} mol")
    print(f"  Theoretical Mass (100% yield): {mass_product_theoretical_g:.2f} g")

# Run the calculation
calculate_reaction_stoichiometry()
