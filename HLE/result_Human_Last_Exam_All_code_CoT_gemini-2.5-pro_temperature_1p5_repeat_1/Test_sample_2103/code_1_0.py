import sys

def find_molecular_formula():
    """
    Determines a molecular formula from high-resolution mass spectrometry data.
    """
    # --- Step 1: Define constants and givens ---
    # High-precision atomic masses (monoisotopic)
    H_MASS = 1.00782503207
    C_MASS = 12.0
    N_MASS = 14.0030740048
    O_MASS = 15.99491461956
    BR_MASS = 78.9183371

    # Given m/z for the protonated ion [M+H]+
    m_plus_h_mass = 1108.70902

    # The isotopic pattern 1:6:15:20:15:6:1 is characteristic of 6 bromine atoms.
    num_br = 6

    # --- Step 2: Calculate the mass of the neutral molecule ---
    neutral_mass = m_plus_h_mass - H_MASS

    # --- Step 3: Calculate the mass of the non-bromine fragment ---
    mass_of_bromines = num_br * BR_MASS
    fragment_target_mass = neutral_mass - mass_of_bromines
    
    # --- Step 4: Systematically search for the fragment's formula ---
    best_match = {
        'formula_str': None,
        'diff': float('inf')
    }
    # Mass tolerance for the search in Daltons (e.g., 3 mDa)
    tolerance = 0.003

    # Nitrogen Rule: The nominal mass of the neutral molecule is floor(1107.69) = 1107 (odd).
    # This implies an odd number of nitrogen atoms in the molecule.
    print("Searching for a molecular formula for a fragment with target mass {:.5f} Da".format(fragment_target_mass))
    print("Constraint: The number of nitrogen atoms must be odd.\n")

    # Search ranges for C, N, O atoms
    for n_count in range(1, 10, 2):  # Odd numbers for N: 1, 3, 5, 7, 9
        for o_count in range(4, 15):
            for c_count in range(25, 45):
                # Calculate the mass contributed by C, N, and O
                cno_mass = (c_count * C_MASS) + (n_count * N_MASS) + (o_count * O_MASS)
                
                # Calculate the remaining mass that must be accounted for by hydrogens
                h_mass_remainder = fragment_target_mass - cno_mass
                
                if h_mass_remainder > 0:
                    # Determine the closest integer number of hydrogens
                    h_count_float = h_mass_remainder / H_MASS
                    h_count_int = round(h_count_float)
                    
                    # Calculate the mass of the proposed CHNO fragment
                    calculated_mass = cno_mass + (h_count_int * H_MASS)
                    mass_diff = abs(calculated_mass - fragment_target_mass)

                    # Check if the match is within the defined tolerance
                    if mass_diff < tolerance:
                        # Sanity check: Degree of Unsaturation (DBE) must be a non-negative integer.
                        # For a stable, closed-shell molecule with an odd number of nitrogens,
                        # the number of hydrogens must also be odd for the DBE to be an integer.
                        if h_count_int > 0 and h_count_int % 2 != 0:
                            dbe = c_count - (h_count_int / 2.0) + (n_count / 2.0) + 1
                            if dbe >= 0 and dbe == int(dbe):
                                # If this is the best match found so far, save it
                                if mass_diff < best_match['diff']:
                                    best_match = {
                                        'c': c_count, 'h': h_count_int, 'n': n_count, 'o': o_count,
                                        'br': num_br,
                                        'mass': calculated_mass,
                                        'diff': mass_diff
                                    }
    
    # --- Step 5: Print the results ---
    if best_match.get('c'):
        c, h, n, o, br = best_match['c'], best_match['h'], best_match['n'], best_match['o'], best_match['br']
        final_formula = f"C{c}H{h}Br{br}N{n}O{o}"
        
        print("Best match found!")
        print("-" * 30)
        print(f"Final Proposed Molecular Formula: {final_formula}")
        print("-" * 30)
        print("Calculation Breakdown:")
        
        c_mass_total = c * C_MASS
        h_mass_total = h * H_MASS
        br_mass_total = br * BR_MASS
        n_mass_total = n * N_MASS
        o_mass_total = o * O_MASS
        total_mass = c_mass_total + h_mass_total + br_mass_total + n_mass_total + o_mass_total
        
        print(f"Mass of {c} Carbon atoms   = {c_mass_total:.5f} Da")
        print(f"Mass of {h} Hydrogen atoms = {h_mass_total:.5f} Da")
        print(f"Mass of {br} Bromine atoms  = {br_mass_total:.5f} Da")
        print(f"Mass of {n} Nitrogen atoms = {n_mass_total:.5f} Da")
        print(f"Mass of {o} Oxygen atoms   = {o_mass_total:.5f} Da")
        print("-" * 30)
        print(f"Total calculated monoisotopic mass (M) = {total_mass:.5f} Da")
        
        # Verify against original data
        m_plus_h_calculated = total_mass + H_MASS
        error_ppm = ((m_plus_h_calculated - m_plus_h_mass) / m_plus_h_mass) * 1e6
        
        print(f"Calculated [M+H]+ mass = {m_plus_h_calculated:.5f} Da")
        print(f"Observed [M+H]+ mass   = {m_plus_h_mass:.5f} Da")
        print(f"Mass Error             = {error_ppm:.2f} ppm")
        
        # Add the final answer tag
        print(f"\n<<<C{c}H{h}Br{br}N{n}O{o}>>>")

    else:
        print("No plausible formula was found within the given constraints and search parameters.")

if __name__ == '__main__':
    find_molecular_formula()