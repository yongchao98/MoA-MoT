import math

def find_molecular_formula():
    """
    Calculates the molecular formula of a neutral species based on high-resolution
    mass spectrometry data of its protonated ion.
    """

    # --- Step 1: Define constants and input data ---
    # Exact masses of the lightest stable isotopes and a hydrogen atom
    MASS = {
        'C': 12.000000,
        'H': 1.007825,
        'N': 14.003074,
        'O': 15.994915,
        'Br': 78.918337,  # ⁷⁹Br
    }
    MASS_H = MASS['H']

    # Input data from the problem
    mz_measured = 1108.70902
    # The isotopic pattern 1:6:15:20:15:6:1 indicates 6 Bromine atoms
    num_br = 6

    print("--- Step-by-step Analysis ---")
    print(f"1. Isotopic pattern suggests {num_br} bromine atoms.")

    # --- Step 2: Calculate the mass of the neutral monoisotopic species ---
    mass_neutral_measured = mz_measured - MASS_H
    print(f"2. Measured m/z of [M+H]⁺ = {mz_measured}")
    print(f"   Calculated neutral mass (M) = {mz_measured} - {MASS_H} = {mass_neutral_measured:.6f} Da")

    # --- Step 3: Calculate the mass of the non-bromine part ---
    mass_br_part = num_br * MASS['Br']
    mass_remainder_target = mass_neutral_measured - mass_br_part
    print(f"3. Mass of {num_br} ⁷⁹Br atoms = {num_br} * {MASS['Br']} = {mass_br_part:.6f} Da")
    print(f"4. Target mass for the remaining C,H,N,O part = {mass_neutral_measured:.6f} - {mass_br_part:.6f} = {mass_remainder_target:.6f} Da")

    # --- Step 4: Systematically search for the formula of the remainder ---
    print("\n5. Searching for the best formula match for the remainder...")
    
    # Define search parameters
    ppm_tolerance = 5.0
    mass_tolerance = (ppm_tolerance / 1_000_000) * mass_remainder_target
    
    # Plausible element ranges for a natural product of this size
    c_range = range(20, 50)
    n_range = range(1, 10)
    o_range = range(1, 15)
    
    best_match = {
        "formula": None,
        "error": float('inf')
    }

    for n_c in c_range:
        for n_n in n_range:
            for n_o in o_range:
                # Calculate the mass filled by C, N, O
                mass_cno = n_c * MASS['C'] + n_n * MASS['N'] + n_o * MASS['O']
                
                # Calculate the remaining mass that must be from H
                mass_h_needed = mass_remainder_target - mass_cno
                
                if mass_h_needed > 0:
                    # Calculate the number of H atoms
                    num_h = mass_h_needed / MASS['H']
                    
                    # Check if the number of H atoms is very close to an integer
                    if abs(num_h - round(num_h)) < 0.02:
                        n_h = int(round(num_h))
                        
                        # Calculate the mass of the candidate remainder formula
                        candidate_mass = n_c * MASS['C'] + n_h * MASS['H'] + n_n * MASS['N'] + n_o * MASS['O']
                        
                        # Check if the mass is within tolerance
                        error = abs(candidate_mass - mass_remainder_target)
                        if error < mass_tolerance and error < best_match["error"]:
                            best_match["formula"] = f"C{n_c}H{n_h}N{n_n}O{n_o}"
                            best_match["error"] = error

    # --- Step 5: Present the final result ---
    if best_match["formula"]:
        # Reconstruct the full formula string
        c_part = best_match["formula"].split('H')[0]
        h_part = 'H' + best_match["formula"].split('H')[1].split('N')[0]
        n_part = 'N' + best_match["formula"].split('N')[1].split('O')[0]
        o_part = 'O' + best_match["formula"].split('O')[1]
        
        final_formula_str = f"{c_part}{h_part}Br{num_br}{n_part}{o_part}"
        
        print(f"\n--- Result ---")
        print(f"The best matching formula for the remainder is: {best_match['formula']}")
        print(f"The molecular formula of the neutral species is: {final_formula_str}")

        # Final verification of the full formula mass
        parts = {'C': int(c_part[1:]), 'H': int(h_part[1:]), 'Br': num_br, 'N': int(n_part[1:]), 'O': int(o_part[1:])}
        
        calculated_mass = sum(parts[el] * MASS[el] for el in ['C', 'H', 'N', 'O']) + parts['Br'] * MASS['Br']
        final_error_ppm = ((calculated_mass - mass_neutral_measured) / mass_neutral_measured) * 1_000_000
        
        print("\nFinal Verification:")
        print(f"Experimental Mass: {mass_neutral_measured:.6f}")
        print(f"Calculated Mass for {final_formula_str}: {calculated_mass:.6f}")
        print(f"Mass Error: {final_error_ppm:.2f} ppm")

    else:
        print("No suitable formula was found within the specified constraints.")

# Execute the function
find_molecular_formula()