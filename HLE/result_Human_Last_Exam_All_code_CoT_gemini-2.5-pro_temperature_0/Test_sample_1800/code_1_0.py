import sys

def calculate_ideal_ni_ce_ratio():
    """
    Calculates the ideal Ni/Ce atomic ratio based on an optimal weight percent
    found in scientific literature for Ni-Ceria catalysts.

    The calculation is based on a 100g sample of the catalyst.
    """
    # --- Step 1: Define constants based on literature and chemical data ---
    # Optimal Ni loading is often cited around 12 wt% for WGS reactions.
    optimal_ni_wt_percent = 12.0

    # Molar masses (g/mol)
    molar_mass_ni = 58.69
    molar_mass_ce = 140.12
    molar_mass_o = 16.00
    molar_mass_ceo2 = molar_mass_ce + 2 * molar_mass_o

    print("--- Finding the Ideal Ni/Ce Atomic Ratio ---")
    print(f"Based on experimental literature, an optimal Ni loading is ~{optimal_ni_wt_percent} wt%.")
    print("This script converts this weight percent to an atomic ratio.\n")

    # --- Step 2: Calculate mass of each component in a 100g sample ---
    mass_ni = optimal_ni_wt_percent
    mass_ceo2 = 100.0 - mass_ni

    print("--- Calculation Steps (for a 100g sample) ---")
    print(f"1. Mass of Ni = {mass_ni:.2f} g")
    print(f"2. Mass of Ceria (CeO2) = {mass_ceo2:.2f} g\n")

    # --- Step 3: Calculate moles of Ni and Ce ---
    moles_ni = mass_ni / molar_mass_ni
    moles_ce = mass_ceo2 / molar_mass_ceo2 # Moles of Ce is the same as moles of CeO2

    print("3. Moles of Ni Calculation:")
    print(f"   Moles Ni = Mass Ni / Molar Mass Ni")
    print(f"   Moles Ni = {mass_ni:.2f} g / {molar_mass_ni} g/mol = {moles_ni:.4f} mol")
    print("\n4. Moles of Ce Calculation:")
    print(f"   Moles Ce = Mass CeO2 / Molar Mass CeO2")
    print(f"   Moles Ce = {mass_ceo2:.2f} g / {molar_mass_ceo2:.2f} g/mol = {moles_ce:.4f} mol")

    # --- Step 4: Calculate the final atomic ratio ---
    atomic_ratio_ni_ce = moles_ni / moles_ce

    print("\n--- Final Atomic Ratio Equation ---")
    print("Ni/Ce Ratio = Moles Ni / Moles Ce")
    # The final print statement shows the numbers in the final equation as requested
    print(f"Ni/Ce Ratio = {moles_ni:.4f} / {moles_ce:.4f}")
    print("\n-----------------------------------------")
    print(f"The ideal Ni/Ce atomic ratio is approximately: {atomic_ratio_ni_ce:.2f}")
    print("-----------------------------------------")
    
    # The final answer is returned in the required format for the system.
    # Using sys.stdout.write to avoid the extra newline from print().
    sys.stdout.write(f"<<<{atomic_ratio_ni_ce:.2f}>>>")


if __name__ == '__main__':
    calculate_ideal_ni_ce_ratio()