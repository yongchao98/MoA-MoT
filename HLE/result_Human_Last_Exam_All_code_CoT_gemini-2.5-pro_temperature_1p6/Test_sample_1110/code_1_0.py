def calculate_reagents(mass_substrate_g, nBuLi_conc_M, nBuLi_eq, borate_eq):
    """
    Calculates the required amount of reagents for the borylation reaction.
    This demonstrates the importance of using precise amounts of n-BuLi.
    """
    # Molecular Weights (g/mol)
    mw_substrate = 317.33  # 2-bromo-4-chloro-1-iodobenzene (C6H3BrClI)
    mw_borate = 103.91     # Trimethyl borate (B(OMe)3)
    density_borate = 0.915 # g/mL

    # --- Calculations ---
    # 1. Moles of the starting material (substrate)
    moles_substrate = mass_substrate_g / mw_substrate

    # 2. Moles of n-BuLi
    moles_nBuLi = moles_substrate * nBuLi_eq

    # 3. Volume of n-BuLi solution needed (in mL)
    # Volume (L) = moles / concentration (mol/L)
    volume_nBuLi_mL = (moles_nBuLi / nBuLi_conc_M) * 1000

    # 4. Amount of trimethyl borate
    moles_borate = moles_substrate * borate_eq
    mass_borate_g = moles_borate * mw_borate
    volume_borate_mL = mass_borate_g / density_borate

    # --- Output ---
    print("Reaction Plan to Avoid Side Products:")
    print("The formation of two boron products suggests a side reaction, likely a second lithiation at the bromine position due to excess n-BuLi.")
    print("To solve this, use a precise, stoichiometric amount of n-BuLi (1.00 eq) after accurately determining its concentration via titration.")
    print("\n--- Example Calculation ---")
    print(f"Starting with {mass_substrate_g:.3f} g of 2-bromo-4-chloro-1-iodobenzene:")
    print(f"This corresponds to {moles_substrate:.4f} moles.")
    print("\nTo add a precise stoichiometric amount:")
    print(f"Target n-BuLi equivalents: {nBuLi_eq:.2f} eq")
    print(f"Required moles of n-BuLi: {moles_nBuLi:.4f} moles")

    # Output the final equation for n-BuLi volume calculation
    print("\nFinal Equation for n-BuLi Volume:")
    print("Vol_nBuLi(mL) = (Mass_Substrate / MW_Substrate) * (Equivalents_nBuLi / Conc_nBuLi) * 1000")
    print(f"Vol_nBuLi(mL) = ({mass_substrate_g:.3f} g / {mw_substrate} g/mol) * ({nBuLi_eq:.2f} eq / {nBuLi_conc_M} M) * 1000 = {volume_nBuLi_mL:.2f} mL")

    print(f"\nAlso add {borate_eq:.1f} eq of trimethyl borate: {volume_borate_mL:.2f} mL ({mass_borate_g:.2f} g)")


if __name__ == '__main__':
    # --- User-defined Parameters ---
    # Mass of starting material in grams
    starting_mass_g = 5.0
    # Concentration of n-BuLi solution in Molarity (mol/L), determined by titration
    nBuLi_concentration = 2.5
    # To solve the problem, we use a more precise amount of n-BuLi
    nBuLi_equivalents = 1.00
    # Equivalents of trimethyl borate
    trimethyl_borate_equivalents = 5.0

    calculate_reagents(starting_mass_g, nBuLi_concentration, nBuLi_equivalents, trimethyl_borate_equivalents)
