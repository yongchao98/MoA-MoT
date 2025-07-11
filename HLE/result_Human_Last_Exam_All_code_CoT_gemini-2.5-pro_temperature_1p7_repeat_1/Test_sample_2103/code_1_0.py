import math

def find_molecular_formula():
    """
    Determines a molecular formula from high-resolution mass spectrometry data
    based on isotopic patterns, exact mass, and chemical principles.
    """
    # --- Step 1: Define constants and input data ---
    # Monoisotopic masses of the most common isotopes
    MASS_H = 1.007825
    MASS_C = 12.000000
    MASS_N = 14.003074
    MASS_O = 15.994915
    MASS_BR = 78.918337  # For 79Br

    # Input data from the problem
    mz_protonated = 1108.70902
    
    # Based on the 1:6:15:20:15:6:1 pattern
    num_br = 6

    # --- Step 2: Calculate the mass of the neutral species ---
    mass_neutral_target = mz_protonated - MASS_H

    # --- Step 3: Calculate the residual mass after subtracting Br6 ---
    mass_br6 = num_br * MASS_BR
    mass_residual = mass_neutral_target - mass_br6
    
    print("--- Analysis Step-by-Step ---")
    print(f"1. Isotopic pattern suggests {num_br} Bromine atoms.")
    print(f"2. Observed m/z for [M+H]+: {mz_protonated:.5f}")
    print(f"3. Calculated neutral monoisotopic mass (M): {mass_neutral_target:.5f} Da")
    print(f"4. Mass of {num_br} Bromine atoms (Br6): {mass_br6:.5f} Da")
    print(f"5. Residual mass for C, H, N, O: {mass_residual:.5f} Da")
    print("\n--- Searching for Formula ---")
    print("Constraints: N is odd, mass error < 5 ppm, Degree of Unsaturation (DoU) is a non-negative integer.\n")

    # --- Step 4: Search for the best combination of C, H, N, O ---
    best_candidate = None
    min_error_ppm = float('inf')

    # Set reasonable search ranges for a complex natural product
    # Nitrogen Rule: Odd mass means an odd number of Nitrogens
    for n_N in range(1, 15, 2):
        for n_O in range(1, 20):
            # Calculate remaining mass for Carbon and Hydrogen
            mass_for_ch = mass_residual - (n_N * MASS_N) - (n_O * MASS_O)
            if mass_for_ch < 0:
                continue

            # Estimate number of carbons and iterate around it
            # Number of hydrogens cannot exceed 2C+2+N for a valid structure
            # A rough guess for n_C is sufficient
            num_c_estimate = int(mass_for_ch / MASS_C)
            
            # Check a small range around the estimate
            for n_C in range(num_c_estimate - 2, num_c_estimate + 3):
                if n_C <= 0:
                    continue
                
                # Calculate required number of hydrogens
                mass_for_h = mass_for_ch - (n_C * MASS_C)
                n_H_float = mass_for_h / MASS_H
                
                # Check if the number of hydrogens is close to an integer
                if abs(n_H_float - round(n_H_float)) < 0.01:
                    n_H = round(n_H_float)
                    if n_H < 0:
                        continue
                        
                    # We have a potential candidate, let's validate it
                    
                    # Calculate total mass
                    calculated_mass = (n_C * MASS_C) + (n_H * MASS_H) + (n_N * MASS_N) + \
                                      (n_O * MASS_O) + mass_br6
                    
                    # Calculate mass error in ppm
                    error_ppm = ((calculated_mass - mass_neutral_target) / mass_neutral_target) * 1e6
                    
                    # Calculate Degrees of Unsaturation (DoU)
                    # DoU = C - H/2 - X/2 + N/2 + 1
                    dou = n_C - (n_H / 2.0) - (num_br / 2.0) + (n_N / 2.0) + 1.0

                    # Check if all conditions are met
                    if abs(error_ppm) < 5 and dou >= 0 and dou == int(dou):
                        # If this is the best hit so far, save it
                        if abs(error_ppm) < abs(min_error_ppm):
                            min_error_ppm = error_ppm
                            best_candidate = {
                                'C': n_C, 'H': n_H, 'N': n_N, 'O': n_O, 'Br': num_br,
                                'mass': calculated_mass, 'error': error_ppm, 'dou': dou
                            }

    # --- Print the final result ---
    if best_candidate:
        c = best_candidate['C']
        h = best_candidate['H']
        n = best_candidate['N']
        o = best_candidate['O']
        br = best_candidate['Br']
        
        print("--- Solution Found ---")
        print(f"The most likely molecular formula is: C{c}H{h}N{n}O{o}Br{br}\n")
        print("Mass calculation for the proposed neutral formula:")
        print(f" {c:2d} * {MASS_C:.6f} (C)  = {c * MASS_C:10.6f} Da")
        print(f" {h:2d} * {MASS_H:.6f} (H)  = {h * MASS_H:10.6f} Da")
        print(f" {n:2d} * {MASS_N:.6f} (N)  = {n * MASS_N:10.6f} Da")
        print(f" {o:2d} * {MASS_O:.6f} (O)  = {o * MASS_O:10.6f} Da")
        print(f" {br:2d} * {MASS_BR:.6f} (Br) = {br * MASS_BR:10.6f} Da")
        print("-------------------------------------------------")
        print(f" Total Calculated Mass   = {best_candidate['mass']:10.6f} Da")
        print(f" Target Mass             = {mass_neutral_target:10.6f} Da")
        print("-------------------------------------------------")
        print(f"Mass Error: {best_candidate['error']:.2f} ppm")
        print(f"Degrees of Unsaturation (DoU): {int(best_candidate['dou'])}")
    else:
        print("No suitable molecular formula was found with the given constraints.")

if __name__ == '__main__':
    find_molecular_formula()
    print("\n<<<C25H25N5O8Br6>>>")