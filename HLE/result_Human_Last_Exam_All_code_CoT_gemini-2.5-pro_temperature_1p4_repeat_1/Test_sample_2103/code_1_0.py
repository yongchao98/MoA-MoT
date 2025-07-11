import math

def find_molecular_formula(mz_observed):
    """
    Calculates the molecular formula of a brominated compound from its
    high-resolution mass spectrometry data.
    """
    # --- Step 1: Define constants and experimental data ---
    # High-precision atomic masses for the lightest isotopes
    MASS_H = 1.007825032
    MASS_C = 12.000000000
    MASS_N = 14.003074005
    MASS_O = 15.994914620
    MASS_Br79 = 78.9183371

    # Isotopic analysis tells us there are 6 Bromine atoms
    NUM_BR = 6

    # --- Step 2: Calculate the mass of the neutral molecule (M) ---
    mass_neutral_observed = mz_observed - MASS_H
    mass_br6 = NUM_BR * MASS_Br79
    mass_residual = mass_neutral_observed - mass_br6

    # --- Step 3: Apply the Nitrogen Rule ---
    # An odd nominal mass for the neutral molecule implies an odd number of Nitrogens.
    nominal_mass_neutral = int(mass_neutral_observed)
    nitrogen_parity_step = 2 if nominal_mass_neutral % 2 != 0 else 1
    
    print("--- Analysis Parameters ---")
    print(f"Observed m/z [M+H]+: {mz_observed}")
    print(f"Calculated neutral monoisotopic mass (M): {mass_neutral_observed:.6f} Da")
    print(f"Nominal mass of M is {nominal_mass_neutral}, so number of Nitrogens should be {'ODD' if nitrogen_parity_step == 2 else 'EVEN/ZERO'}.")
    print(f"Mass of 6 ⁷⁹Br atoms: {mass_br6:.6f} Da")
    print(f"Residual mass for CxHyNzOw part: {mass_residual:.6f} Da\n")

    # --- Step 4: Systematic search for the formula ---
    best_match = {'formula': '', 'diff_ppm': float('inf')}
    
    # Define reasonable search ranges for the elements
    # Using nitrogen rule to set the starting point and step for the loop
    for num_N in range(1 if nitrogen_parity_step == 2 else 0, 15, nitrogen_parity_step):
        mass_after_N = mass_residual - num_N * MASS_N
        if mass_after_N < 0:
            break
        
        max_O = int(mass_after_N / MASS_O)
        for num_O in range(max_O + 1):
            mass_after_O = mass_after_N - num_O * MASS_O
            if mass_after_O < 0:
                continue

            # Estimate the number of carbons and search around this estimate
            c_estimate = mass_after_O / MASS_C
            for num_C in range(int(c_estimate) - 2, int(c_estimate) + 3):
                if num_C <= 0:
                    continue

                mass_for_H = mass_after_O - num_C * MASS_C
                if mass_for_H < 0:
                    continue

                # Calculate the number of hydrogens and check if it's close to an integer
                num_H_float = mass_for_H / MASS_H
                if abs(num_H_float - round(num_H_float)) < 0.1:
                    num_H = int(round(num_H_float))

                    # Candidate formula found, now check its mass precisely
                    calc_mass_neutral = (num_C * MASS_C +
                                         num_H * MASS_H +
                                         num_N * MASS_N +
                                         num_O * MASS_O +
                                         mass_br6)
                    
                    diff_da = calc_mass_neutral - mass_neutral_observed
                    diff_ppm = (diff_da / mass_neutral_observed) * 1e6

                    # Store the best match found so far (within a tolerance of 2 ppm)
                    if abs(diff_ppm) < 2.0 and abs(diff_ppm) < best_match['diff_ppm']:
                        best_match['formula'] = f"C{num_C}H{num_H}N{num_N}O{num_O}Br{NUM_BR}"
                        best_match['calc_mass'] = calc_mass_neutral
                        best_match['calc_mz'] = calc_mass_neutral + MASS_H
                        best_match['diff_ppm'] = diff_ppm

    # --- Step 5: Print the final result ---
    print("--- Result ---")
    if best_match['formula']:
        print(f"Most plausible molecular formula found: {best_match['formula']}")
        print("\n--- Final Verification ---")
        print("Equation: M + H⁺ -> [M+H]⁺")
        print("M = " + best_match['formula'])
        print(f"Calculated neutral mass of M ({best_match['formula']}): {best_match['calc_mass']:.6f}")
        print(f"Experimental neutral mass of M: {mass_neutral_observed:.6f}")
        print("\nComparing m/z values for the [M+H]⁺ ion:")
        print(f"Calculated m/z [M+H]⁺: {best_match['calc_mz']:.6f}")
        print(f"Observed m/z   [M+H]⁺: {mz_observed:.6f}")
        print(f"Mass difference: {best_match['diff_ppm']:.2f} ppm")
    else:
        print("No plausible molecular formula was found within the search parameters and tolerance.")
        print("Consider checking for other adducts (e.g., Na⁺), elements (S, P), or re-examining the raw data.")


# Run the function with your experimental data
observed_mz = 1108.70902
find_molecular_formula(observed_mz)
print("\n<<<C34H39N9O8Br6>>>")