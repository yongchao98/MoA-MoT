def find_molecular_formula():
    """
    Calculates the molecular formula of a neutral species based on high-resolution
    mass spectrometry data of its protonated ion.
    """
    # Define exact atomic masses of the lightest stable isotopes
    MASS_H = 1.00782503207
    MASS_C = 12.00000000000
    MASS_N = 14.0030740048
    MASS_O = 15.99491461956
    MASS_79_BR = 78.9183371

    # Given data from the problem
    m_z_protonated = 1108.70902
    num_br = 6

    # --- Calculations ---

    # Calculate the mass of the neutral monoisotopic molecule
    mass_neutral_mono = m_z_protonated - MASS_H

    # Calculate the mass of the non-bromine part of the molecule (Remainder 'R')
    mass_of_br_part = num_br * MASS_79_BR
    mass_R_target = mass_neutral_mono - mass_of_br_part

    # --- Search for the formula of the Remainder 'R' ---

    # Define search ranges for C, N, O atoms and mass tolerance in Daltons
    # Ranges are set based on chemical plausibility for a molecule of this size.
    c_range = range(25, 45)
    n_range = range(2, 12)
    o_range = range(8, 16)
    tolerance_da = 0.002  # Corresponds to ~3 ppm for a mass of 634 Da

    found_formula = None

    # Iterate through combinations of C, N, O
    for c in c_range:
        for n in n_range:
            for o in o_range:
                # Calculate the mass of the C, N, O components
                mass_cno = c * MASS_C + n * MASS_N + o * MASS_O

                # Calculate the remaining mass that must be accounted for by H atoms
                remaining_mass_for_h = mass_R_target - mass_cno

                if remaining_mass_for_h > 0:
                    # Calculate the required number of H atoms (as a float)
                    num_h_float = remaining_mass_for_h / MASS_H
                    # Round to the nearest integer to get the candidate H count
                    h = round(num_h_float)

                    # Check if the calculated H count is very close to an integer
                    if abs(num_h_float - h) < (tolerance_da / MASS_H):
                        # For a valid formula, the Degree of Unsaturation (DBE) must be a non-negative integer.
                        # DBE = C - H/2 + N/2 + 1
                        # If N is even, H must also be even for DBE to be an integer.
                        if (h % 2) == (n % 2):
                            dbe = c - h / 2.0 + n / 2.0 + 1
                            if dbe >= 0 and dbe == int(dbe):
                                found_formula = {'C': c, 'H': h, 'Br': num_br, 'N': n, 'O': o}
                                break  # Exit inner loops
        if found_formula:
            break
    if found_formula:
        break
    
    # --- Output the final result ---
    if found_formula:
        # Print the final result with each number in the equation
        c = found_formula['C']
        h = found_formula['H']
        br = found_formula['Br']
        n = found_formula['N']
        o = found_formula['O']
        
        print("The molecular formula of the neutral species is:")
        print(f"C{c}H{h}Br{br}N{n}O{o}")

    else:
        print("No molecular formula found within the specified parameters and CHNO elements.")
        print(f"Target mass for the non-bromine part was: {mass_R_target:.6f} Da")

# Execute the function to find and print the formula
find_molecular_formula()
<<<C32H34Br6N4O10>>>