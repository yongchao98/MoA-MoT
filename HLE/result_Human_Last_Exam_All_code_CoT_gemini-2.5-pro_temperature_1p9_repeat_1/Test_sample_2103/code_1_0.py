import math

def find_molecular_formula():
    """
    Calculates the molecular formula of a neutral species based on high-resolution mass spectrometry data.
    """
    # --- Step 1: Define constants and initial data ---
    mz_protonated_ion = 1108.70902
    num_bromine = 6
    
    # Monoisotopic masses of the most common stable isotopes
    MASS_H = 1.007825032
    MASS_C = 12.000000000
    MASS_N = 14.003074004
    MASS_O = 15.994914620
    MASS_BR79 = 78.9183371

    print("--- Step-by-step analysis ---")
    
    # --- Step 2: Calculate the mass of the neutral molecule (M) ---
    mass_neutral_M = mz_protonated_ion - MASS_H
    print(f"1. Isotopic pattern (1:6:15:20:15:6:1) indicates the presence of {num_bromine} bromine atoms.")
    print(f"2. Mass of the protonated ion [M+H]+ is {mz_protonated_ion}.")
    print(f"3. Calculated monoisotopic mass of the neutral molecule (M) is {mass_neutral_M:.6f} Da.")

    # --- Step 3: Calculate the mass of the non-bromine part ---
    mass_of_6_br = num_bromine * MASS_BR79
    target_mass = mass_neutral_M - mass_of_6_br
    print(f"4. Mass of {num_bromine} ⁷⁹Br atoms is {mass_of_6_br:.6f} Da.")
    print(f"5. Target mass for the remaining CxHyNzOw formula is {target_mass:.6f} Da.")
    print("\n--- Searching for the molecular formula ---")

    # --- Step 4: Systematically search for the formula ---
    # Search ranges for the number of atoms
    max_C = 50
    max_O = 20
    # Nitrogen Rule: Mass(M) is odd (1107), so N must be odd.
    max_N = 15 
    
    # Mass tolerance for a match (e.g., 5 ppm)
    tolerance = 5e-6 * target_mass 
    
    best_match = None
    smallest_error = float('inf')

    # Start search. Nitrogen is odd.
    for n_N in range(1, max_N, 2):
        for n_C in range(1, max_C):
            for n_O in range(max_O + 1):
                # For a stable molecule, 2*DBE = 2*C - H + N + 2 must be an even non-negative integer.
                # Since N is odd, H must also be odd for -H+N to be even.
                
                # Calculate the mass required for H atoms
                current_mass = (n_C * MASS_C) + (n_N * MASS_N) + (n_O * MASS_O)
                rem_mass_for_H = target_mass - current_mass
                
                # Estimate the number of H atoms
                if rem_mass_for_H < 0:
                    continue
                n_H_float = rem_mass_for_H / MASS_H
                n_H = round(n_H_float)
                
                # H must be a positive odd integer
                if n_H <= 0 or n_H % 2 == 0:
                    continue

                # Check how close the calculated n_H is to an integer. If it's not very close, it's not a good fit.
                if abs(n_H - n_H_float) > 0.01:
                    continue

                # Calculate the full mass and error for the candidate formula
                calculated_mass = current_mass + n_H * MASS_H
                error = abs(calculated_mass - target_mass)

                # Check for plausibility (DBE must be a non-negative integer)
                dbe = n_C - n_H / 2 + n_N / 2 + 1
                if dbe >= 0 and dbe == int(dbe):
                    if error < smallest_error:
                        smallest_error = error
                        best_match = {'C': n_C, 'H': n_H, 'N': n_N, 'O': n_O, 'Br': num_bromine}
    
    if best_match and smallest_error < tolerance:
        formula = best_match
        print(f"\nFound a plausible formula: C{formula['C']}H{formula['H']}N{formula['N']}O{formula['O']}Br{formula['Br']}")

        # --- Final verification and output ---
        print("\n--- Verifying the final molecular formula ---")
        mass_C_total = formula['C'] * MASS_C
        mass_H_total = formula['H'] * MASS_H
        mass_N_total = formula['N'] * MASS_N
        mass_O_total = formula['O'] * MASS_O
        mass_Br_total = formula['Br'] * MASS_BR79
        
        total_calc_mass = mass_C_total + mass_H_total + mass_N_total + mass_O_total + mass_Br_total
        
        print("Final calculation breakdown:")
        print(f"C{formula['C']}: {formula['C']:>2} * {MASS_C:12.7f} = {mass_C_total:9.6f} Da")
        print(f"H{formula['H']}: {formula['H']:>2} * {MASS_H:12.7f} = {mass_H_total:9.6f} Da")
        print(f"N{formula['N']}: {formula['N']:>2} * {MASS_N:12.7f} = {mass_N_total:9.6f} Da")
        print(f"O{formula['O']}: {formula['O']:>2} * {MASS_O:12.7f} = {mass_O_total:9.6f} Da")
        print(f"Br{formula['Br']}: {formula['Br']:>2} * {MASS_BR79:12.7f} = {mass_Br_total:9.6f} Da")
        print("-------------------------------------------------")
        print(f"Total calculated mass:      {total_calc_mass:15.6f} Da")
        print(f"Experimental neutral mass:  {mass_neutral_M:15.6f} Da")
        print(f"Error:                      {abs(total_calc_mass - mass_neutral_M):15.6f} Da")
    else:
        print("Could not find a matching formula within the specified parameters.")

find_molecular_formula()