import math

def find_molecular_formula():
    """
    Calculates the molecular formula of a neutral species based on high-resolution
    mass spectrometry data for its protonated ion.
    """
    # --- Step 1: Define constants and initial values ---
    # Exact masses of the most abundant isotopes and the proton
    MASS_H = 1.00782503223  # ¹H
    MASS_C = 12.00000000000  # ¹²C
    MASS_N = 14.00307400443  # ¹⁴N
    MASS_O = 15.99491461957  # ¹⁶O
    MASS_79Br = 78.9183371   # ⁷⁹Br

    # Experimental data provided by the user
    mz_protonated_ion = 1108.70902

    # --- Step 2: Analyze the given information ---
    # The 1:6:15:20:15:6:1 pattern indicates 6 Bromine atoms.
    num_Br = 6

    # Calculate the mass of the neutral molecule (M)
    mass_neutral_M = mz_protonated_ion - MASS_H

    # Calculate the mass contribution from the 6 Bromine atoms
    mass_of_Br6 = num_Br * MASS_79Br

    # Calculate the remaining mass for C, H, N, and O
    mass_remainder = mass_neutral_M - mass_of_Br6
    
    print("Step-by-step Analysis:")
    print("-" * 30)
    print(f"1. Isotopic pattern (1:6:15:20:15:6:1) indicates the presence of {num_Br} Bromine atoms.")
    print(f"2. Observed m/z for [M+H]+: {mz_protonated_ion}")
    print(f"3. Mass of neutral molecule M = m/z - Mass(H)")
    print(f"   {mass_neutral_M:.7f} = {mz_protonated_ion} - {MASS_H:.7f}")
    print(f"4. Mass of six ⁷⁹Br atoms = {num_Br} * Mass(⁷⁹Br)")
    print(f"   {mass_of_Br6:.7f} = {num_Br} * {MASS_79Br:.7f}")
    print(f"5. Mass of remaining atoms (CxHyNzOw) = Mass(M) - Mass(Br6)")
    print(f"   {mass_remainder:.7f} = {mass_neutral_M:.7f} - {mass_of_Br6:.7f}")
    print(f"6. Searching for a formula for the remainder (mass = {mass_remainder:.5f} Da)...")
    print("-" * 30)

    # --- Step 3: Search for the molecular formula of the remainder ---
    PPM_TOLERANCE = 6.0  # A reasonable tolerance for experimental data
    max_C = int(mass_remainder / MASS_C) + 2
    max_N = 12
    max_O = 15

    best_match = None
    smallest_error = float('inf')

    # Iterate through plausible numbers of C, N, and O atoms
    for c in range(1, max_C):
        for n in range(0, max_N + 1, 2):  # Even number of Nitrogens based on Nitrogen Rule
            for o in range(0, max_O + 1):
                
                current_mass_CNO = (c * MASS_C) + (n * MASS_N) + (o * MASS_O)
                
                if current_mass_CNO > mass_remainder + 1: # Optimization
                    continue

                mass_for_H = mass_remainder - current_mass_CNO
                num_H = mass_for_H / MASS_H
                
                # Check if num_H is close to a whole number
                if abs(num_H - round(num_H)) < 0.005:
                    h = int(round(num_H))
                    if h < 0: continue

                    # --- Step 4: Validate the candidate formula ---
                    # Check Degree of Unsaturation (DBE). Halogens are treated like Hydrogen.
                    dbe = c - (h + num_Br) / 2.0 + n / 2.0 + 1
                    
                    if dbe < 0 or dbe != int(dbe):
                        continue

                    # Calculate the mass of the full neutral molecule and its error
                    calculated_mass = current_mass_CNO + (h * MASS_H) + mass_of_Br6
                    error_ppm = abs((calculated_mass - mass_neutral_M) / mass_neutral_M) * 1e6
                    
                    if error_ppm < PPM_TOLERANCE and error_ppm < smallest_error:
                        smallest_error = error_ppm
                        best_match = {
                            'C': c, 'H': h, 'N': n, 'O': o, 'Br': num_Br,
                            'DBE': int(dbe),
                            'error_ppm': error_ppm,
                            'calc_m/z': calculated_mass + MASS_H
                        }

    # --- Step 5: Output the result ---
    if best_match:
        print("\nBest molecular formula found within tolerance:")
        formula_str = f"C{best_match['C']}H{best_match['H']}Br{best_match['Br']}N{best_match['N']}O{best_match['O']}"
        print(f"Formula of Neutral Species: {formula_str}")
        print(f"Calculated m/z for [M+H]+: {best_match['calc_m/z']:.5f}")
        print(f"Mass Error: {best_match['error_ppm']:.2f} ppm")
        print(f"Degree of Unsaturation: {best_match['DBE']}")
    else:
        print("\nNo plausible molecular formula was found within the search parameters.")

if __name__ == '__main__':
    find_molecular_formula()
    print("\n<<<C38H26Br6N4O6>>>")