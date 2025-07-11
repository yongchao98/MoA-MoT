import math

def find_molecular_formula():
    """
    Calculates the molecular formula of a neutral species based on high-resolution mass spectrometry data.
    """
    # --- Part 1: Define constants and input data ---

    # Exact monoisotopic masses of elements (Da)
    MASS = {
        'H': 1.007825032,
        'C': 12.000000000,
        'N': 14.003074004,
        'O': 15.994914620,
        'Br79': 78.9183371
    }
    
    # Given data from the problem
    mz_m_plus_h = 1108.70902
    num_br = 6
    ppm_tolerance = 5.0 # A realistic tolerance for high-resolution mass spec

    # --- Part 2: Calculate target masses ---

    # Calculate the monoisotopic mass of the neutral species (M)
    mass_neutral = mz_m_plus_h - MASS['H']

    # Calculate the mass contribution of the six 79Br atoms
    mass_br_part = num_br * MASS['Br79']
    
    # Calculate the residual mass to be explained by C, H, N, and O
    mass_residual_target = mass_neutral - mass_br_part
    
    print("--- Analysis Step-by-Step ---")
    print(f"1. Isotopic pattern suggests {num_br} Bromine atoms.")
    print(f"2. Observed m/z for [M+H]+: {mz_m_plus_h:.5f}")
    print(f"3. Calculated neutral monoisotopic mass (M): {mass_neutral:.5f} Da")
    print(f"4. Mass of Br_6 part: {mass_br_part:.5f} Da")
    print(f"5. Target mass for C,H,N,O part: {mass_residual_target:.5f} Da\n")
    print("--- Searching for C,H,N,O Composition ---")
    print("Constraints: N must be odd, H must be odd, DoU must be a non-negative integer.\n")

    # --- Part 3: Brute-force search for the C,H,N,O composition ---

    # Set reasonable search ranges for natural products of this size
    # Iterate heavier atoms first to prune the search space more quickly
    for num_o in range(1, 25):
        for num_c in range(10, 60):
            # Nitrogen Rule implies an odd number of N atoms
            for num_n in range(1, 15, 2):
                
                # Calculate the mass of the C, N, and O atoms in the current guess
                mass_c_n_o = num_c * MASS['C'] + num_n * MASS['N'] + num_o * MASS['O']
                
                # Calculate the remaining mass that must be accounted for by Hydrogen
                mass_for_h = mass_residual_target - mass_c_n_o

                # If the mass of C,N,O alone is already too high, skip to the next iteration
                if mass_for_h < 0:
                    continue
                
                # Estimate the number of H atoms required
                num_h = round(mass_for_h / MASS['H'])

                # Apply Parity Rule: H must be odd since N is odd. Skip if num_h is non-positive or even.
                if num_h <= 0 or num_h % 2 == 0:
                    continue
                
                # Calculate the exact mass of the C,H,N,O candidate formula
                candidate_residual_mass = mass_c_n_o + num_h * MASS['H']
                
                # Check if the candidate's mass is within the ppm tolerance of the target
                ppm_error = ((candidate_residual_mass - mass_residual_target) / mass_residual_target) * 1e6
                
                if abs(ppm_error) < ppm_tolerance:
                    # Potential hit! Now, check for chemical feasibility using Degrees of Unsaturation (DoU).
                    # DoU = C - H/2 - X/2 + N/2 + 1  (where X is number of halogens)
                    dou = num_c - (num_h / 2.0) - (num_br / 2.0) + (num_n / 2.0) + 1
                    
                    # DoU should be a non-negative integer for a stable molecule.
                    if dou >= 0 and abs(dou - round(dou)) < 0.001:
                        dou_int = int(round(dou))
                        
                        # --- Part 4: Found a valid solution. Print the results. ---
                        print(">>> Found a matching formula! <<<\n")
                        final_formula = f"C{num_c}H{num_h}N{num_n}O{num_o}Br{num_br}"
                        print(f"Proposed Molecular Formula: {final_formula}")
                        
                        calculated_total_mass = candidate_residual_mass + mass_br_part
                        
                        print("\n--- Verification of Mass ---")
                        print(f"  {num_c:<2} * {MASS['C']:<12.7f} (C)  = {num_c * MASS['C']:>10.5f}")
                        print(f"  {num_h:<2} * {MASS['H']:<12.7f} (H)  = {num_h * MASS['H']:>10.5f}")
                        print(f"  {num_n:<2} * {MASS['N']:<12.7f} (N)  = {num_n * MASS['N']:>10.5f}")
                        print(f"  {num_o:<2} * {MASS['O']:<12.7f} (O)  = {num_o * MASS['O']:>10.5f}")
                        print(f"  {num_br:<2} * {MASS['Br79']:<12.7f} (Br) = {num_br * MASS['Br79']:>10.5f}")
                        print("  -------------------------------------------")
                        print(f"  Calculated Neutral Mass     = {calculated_total_mass:>10.5f}")
                        print(f"  Experimental Neutral Mass   = {mass_neutral:>10.5f}\n")
                        
                        print("--- Verification of other properties ---")
                        print(f"Mass Error: {ppm_error:.2f} ppm (within {ppm_tolerance} ppm tolerance)")
                        print(f"Degrees of Unsaturation (DoU): {dou_int}")
                        
                        # Return the found formula
                        return final_formula
                        
    return "No valid formula found within the given parameters."

# Run the function and capture the final answer
final_answer = find_molecular_formula()
print(f"\n<<<C{final_answer[1:3]}H{final_answer[4:6]}N{final_answer[7]}O{final_answer[9]}Br{final_answer[12]}>>>")
