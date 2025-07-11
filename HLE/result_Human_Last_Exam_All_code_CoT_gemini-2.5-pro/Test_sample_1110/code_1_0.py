import sys

def calculate_borylation_reagents():
    """
    Calculates reagent quantities for the borylation of 2-bromo-4-chloro-1-iodobenzene
    and explains the importance of precise stoichiometry.
    """
    # --- Constants ---
    # Molecular Weight of 2-bromo-4-chloro-1-iodobenzene (C6H3BrClI)
    MW_SM = (12.01 * 6) + (1.01 * 3) + 79.90 + 35.45 + 126.90  # g/mol
    # Molecular Weight of Trimethyl borate B(OMe)3
    MW_BORATE = 10.81 + (12.01 * 3) + (1.01 * 9) + (16.00 * 3) # g/mol
    DENSITY_BORATE = 0.932  # g/mL

    # --- User Inputs (with default values for demonstration) ---
    mass_sm_g = 10.0  # grams of starting material
    molarity_nBuLi = 2.5  # Molarity of n-BuLi solution (mol/L)

    print("### Stoichiometry Calculation for Borylation Reaction ###\n")
    print(f"This script calculates the reagents needed for the reaction:")
    print("2-bromo-4-chloro-1-iodobenzene + n-BuLi -> (2-bromo-4-chlorophenyl)lithium")
    print("(2-bromo-4-chlorophenyl)lithium + B(OMe)3 -> Target Product\n")
    print("--- Problem Analysis ---")
    print("The observation of two Boron (B) NMR signals means a side-product was formed.")
    print("This is likely a diaryl species (Ar2-BOH), formed when the aryllithium (Ar-Li)")
    print("intermediate attacks the desired boronate product (Ar-B(OMe)2).")
    print("This happens when there is an excess of Ar-Li.")
    print("Solution: Use a precise amount of n-BuLi to avoid generating excess Ar-Li.\n")


    # --- Calculations ---
    moles_sm = mass_sm_g / MW_SM

    # Case 1: Standard protocol with slight excess
    eq_nBuLi_standard = 1.05
    moles_nBuLi_standard = moles_sm * eq_nBuLi_standard
    volume_nBuLi_standard_mL = (moles_nBuLi_standard / molarity_nBuLi) * 1000

    # Case 2: Precise protocol to avoid side products
    eq_nBuLi_precise = 1.00
    moles_nBuLi_precise = moles_sm * eq_nBuLi_precise
    volume_nBuLi_precise_mL = (moles_nBuLi_precise / molarity_nBuLi) * 1000

    # Borate calculation (remains the same)
    eq_borate = 5.0
    moles_borate = moles_sm * eq_borate
    mass_borate_g = moles_borate * MW_BORATE
    volume_borate_mL = mass_borate_g / DENSITY_BORATE

    # --- Output ---
    print("--- Reagent Calculation Example ---")
    print(f"Starting with {mass_sm_g:.2f} g of 2-bromo-4-chloro-1-iodobenzene (MW: {MW_SM:.2f} g/mol):")
    print(f"Moles of Starting Material = {mass_sm_g:.2f} g / {MW_SM:.2f} g/mol = {moles_sm:.4f} mol\n")

    print(f"To add Trimethyl Borate ({eq_borate:.1f} eq):")
    print(f"Volume = {volume_borate_mL:.2f} mL")

    print("--- Comparing n-BuLi Amounts ---")
    print(f"Original Method (1.05 eq n-BuLi):")
    print(f"Volume of {molarity_nBuLi:.2f} M n-BuLi = {volume_nBuLi_standard_mL:.2f} mL")
    print("This slight excess can lead to the diaryl side-product, causing the 2nd NMR signal.\n")

    print(f"Proposed Solution (1.00 eq n-BuLi):")
    print(f"Volume of {molarity_nBuLi:.2f} M n-BuLi = {volume_nBuLi_precise_mL:.2f} mL")
    print("Using a precise amount minimizes excess reactive intermediate, leading to a cleaner reaction.")

    # Final reaction equation summary
    print("\n--- Final Recommended Stoichiometry ---")
    print(f"1.00 eq (2-bromo-4-chloro-1-iodobenzene) : {moles_sm:.4f} mol")
    print(f"1.00 eq (n-BuLi)                       : {moles_nBuLi_precise:.4f} mol")
    print(f"5.00 eq (Trimethyl borate)             : {moles_borate:.4f} mol")


if __name__ == '__main__':
    calculate_borylation_reagents()