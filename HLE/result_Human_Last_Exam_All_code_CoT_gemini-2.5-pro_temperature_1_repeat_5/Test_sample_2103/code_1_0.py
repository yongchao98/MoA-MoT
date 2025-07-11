import math

def find_molecular_formula():
    """
    Calculates the molecular formula of a neutral species based on high-resolution
    mass spectrometry data and isotopic patterns.
    """
    # --- Step 0: Define constants and input data ---
    # Exact masses of the most abundant isotopes
    H_MASS = 1.007825
    C_MASS = 12.000000
    N_MASS = 14.003074
    O_MASS = 15.994915
    BR_MASS = 78.918337  # 79Br

    # Input data from the problem
    mz_protonated = 1108.70902
    # The isotopic pattern 1:6:15:20:15:6:1 indicates 6 bromine atoms
    num_br = 6
    # Mass tolerance in Daltons for finding a match
    mass_tolerance = 0.005

    print("Determining the molecular formula based on the provided MS data...")
    print("-" * 60)

    # --- Step 1: Calculate the mass of the neutral monoisotopic molecule (M) ---
    mass_M = mz_protonated - H_MASS
    print(f"Step 1: Calculate the mass of the neutral molecule.")
    print(f"Mass = m/z_protonated - Mass_H = {mz_protonated} - {H_MASS} = {mass_M:.5f} Da")

    # --- Step 2: Calculate the mass of the non-bromine part (CxHyNzOw) ---
    mass_6Br = num_br * BR_MASS
    mass_remainder = mass_M - mass_6Br
    print(f"\nStep 2: Subtract the mass of the 6 bromine atoms.")
    print(f"Mass_remainder = Mass_neutral - Mass_6Br = {mass_M:.5f} - (6 * {BR_MASS}) = {mass_remainder:.5f} Da")

    print(f"\nStep 3: Computationally search for the elemental composition of the remainder.")
    print("This search is constrained by the Nitrogen Rule (odd N atoms) and mass tolerance.")
    
    # --- Step 3: Iterate through plausible numbers of C, H, N, O ---
    # Set reasonable search ranges for a molecule of this size
    max_N = 15  # Search odd numbers up to 15
    max_O = 20
    max_C = int(mass_remainder / C_MASS) + 2 # Upper limit for Carbon atoms

    found_formula = None

    for n_count in range(1, max_N, 2):  # Nitrogen Rule: odd numbers only
        mass_after_N = mass_remainder - (n_count * N_MASS)
        if mass_after_N < 0: continue

        for o_count in range(max_O + 1):
            mass_after_O = mass_after_N - (o_count * O_MASS)
            if mass_after_O < 0: continue
            
            # From the remaining mass, calculate the required number of carbons
            # mass_after_O = c_count * C_MASS + h_count * H_MASS
            # Since H_MASS is slightly > 1 and C_MASS is 12, h_count ~ mass_after_O % 12
            # A more robust way is to iterate C and calculate H.
            
            c_count = int(round(mass_after_O / C_MASS))
            
            for c_candidate in range(c_count - 2, c_count + 3):
                if c_candidate <= 0: continue
                
                mass_for_H = mass_after_O - (c_candidate * C_MASS)
                if mass_for_H < 0: continue
                
                h_count_float = mass_for_H / H_MASS
                h_count = int(round(h_count_float))

                # Check if the calculated H count is a positive integer within tolerance
                if h_count > 0 and abs(h_count_float - h_count) < (mass_tolerance / H_MASS):
                    
                    # Verify that the total calculated mass is within tolerance
                    calculated_mass = (c_candidate * C_MASS + h_count * H_MASS +
                                       n_count * N_MASS + o_count * O_MASS + mass_6Br)
                    
                    if abs(calculated_mass - mass_M) < mass_tolerance:
                        found_formula = {
                            'C': c_candidate,
                            'H': h_count,
                            'N': n_count,
                            'O': o_count,
                            'Br': num_br,
                            'mass': calculated_mass
                        }
                        break
            if found_formula: break
        if found_formula: break

    # --- Step 4: Output the result ---
    print("\nStep 4: Combine the parts to form the final molecular formula.")
    if found_formula:
        c = found_formula['C']
        h = found_formula['H']
        n = found_formula['N']
        o = found_formula['O']
        br = found_formula['Br']
        mass = found_formula['mass']
        
        print("\nA potential molecular formula has been found:")
        print(f"  - Carbon atoms (C):   {c}")
        print(f"  - Hydrogen atoms (H): {h}")
        print(f"  - Nitrogen atoms (N): {n}")
        print(f"  - Oxygen atoms (O):   {o}")
        print(f"  - Bromine atoms (Br): {br}")

        print("\nVerification:")
        print(f"Calculated Mass for C{c}H{h}N{n}O{o}Br{br}: {mass:.5f} Da")
        print(f"Observed Mass of Neutral Species:     {mass_M:.5f} Da")
        print(f"Mass Difference:                      {abs(mass - mass_M):.5f} Da")
        
        final_formula_str = f"C{c}H{h}N{n}O{o}Br{br}"
        print(f"\nThe molecular formula of the neutral species is: {final_formula_str}")
        return final_formula_str
    else:
        print("\nCould not find a matching molecular formula within the given constraints.")
        return None

# Execute the function and capture the result for the final answer format
final_formula = find_molecular_formula()
if final_formula:
    print(f"\n<<<{final_formula}>>>")
