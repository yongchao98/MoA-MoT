import math

def find_molecular_formula():
    """
    Determines the molecular formula of a neutral species based on high-resolution
    mass spectrometry data.
    """
    # --- Step 0: Define constants and input data ---
    # Precise monoisotopic masses (in atomic mass units, u)
    MASS_C = 12.000000
    MASS_H = 1.007825
    MASS_N = 14.003074
    MASS_O = 15.994915
    MASS_79Br = 78.918337
    MASS_PROTON = 1.007276  # Mass of H+

    # Given experimental data
    mz_observed = 1108.70902
    
    print("### Analysis of the Molecular Formula ###\n")

    # --- Step 1: Interpret the isotopic pattern ---
    num_br = 6
    print(f"Step 1: Isotopic Pattern Analysis")
    print(f"The 1:6:15:20:15:6:1 pattern indicates the presence of {num_br} Bromine atoms.\n")

    # --- Step 2: Calculate the neutral monoisotopic mass ---
    mass_neutral_total = mz_observed - MASS_PROTON
    print(f"Step 2: Calculate Neutral Monoisotopic Mass")
    print(f"Mass of [M+H]+ ion: {mz_observed:.5f} u")
    print(f"Mass of neutral molecule (M) = {mz_observed:.5f} - {MASS_PROTON:.5f} = {mass_neutral_total:.5f} u\n")

    # --- Step 3: Determine the mass of the non-bromine portion ---
    mass_br_part = num_br * MASS_79Br
    target_mass_remainder = mass_neutral_total - mass_br_part
    print(f"Step 3: Calculate Mass of the Remainder (C, H, N, O)")
    print(f"Mass of 6 Bromine atoms (6 * {MASS_79Br:.5f}) = {mass_br_part:.5f} u")
    print(f"Mass of remainder = {mass_neutral_total:.5f} - {mass_br_part:.5f} = {target_mass_remainder:.5f} u\n")

    # --- Step 4: Search for the remaining formula ---
    print(f"Step 4: Searching for the C, H, N, O Formula")
    print("Applying Nitrogen Rule (odd number of Nitrogens for odd nominal mass) and searching...\n")
    
    # Search parameters
    tolerance = 0.003  # Da
    found_formula = None

    # Search loops (ranges based on chemical intuition for a molecule of this size)
    for n_o in range(5, 12):
        for n_n in range(1, 10, 2):  # Odd numbers of Nitrogen
            for n_c in range(30, 40):
                # Calculate mass of C, N, O part
                mass_cno = n_c * MASS_C + n_n * MASS_N + n_o * MASS_O
                
                # Calculate remaining mass for H
                rem_mass_for_h = target_mass_remainder - mass_cno
                
                if rem_mass_for_h > 0:
                    n_h = round(rem_mass_for_h / MASS_H)
                    
                    # Calculate total mass of the C,H,N,O part
                    calculated_mass = mass_cno + n_h * MASS_H
                    
                    # Check if within tolerance
                    if abs(calculated_mass - target_mass_remainder) < tolerance:
                        # Check for valid Degree of Unsaturation (must be >= 0 and half-integer)
                        dbe = n_c - n_h / 2.0 + n_n / 2.0 + 1
                        if dbe >= 0 and abs(dbe - round(dbe, 1)) < 0.01:
                            found_formula = {'C': n_c, 'H': n_h, 'Br': num_br, 'N': n_n, 'O': n_o}
                            break
            if found_formula: break
        if found_formula: break

    # --- Step 5: Final Verification and Output ---
    if found_formula:
        c, h, br, n, o = found_formula['C'], found_formula['H'], found_formula['Br'], found_formula['N'], found_formula['O']
        
        final_formula_str = f"C{c}H{h}Br{br}N{n}O{o}"
        print(f"Found a plausible formula: {final_formula_str}\n")
        
        print("### Final Verification ###")
        mass_C_part = c * MASS_C
        mass_H_part = h * MASS_H
        mass_N_part = n * MASS_N
        mass_O_part = o * MASS_O
        
        recalculated_mass = mass_C_part + mass_H_part + mass_N_part + mass_O_part + mass_br_part
        recalculated_mz = recalculated_mass + MASS_PROTON
        error_ppm = ((recalculated_mz - mz_observed) / mz_observed) * 1e6

        print(f"Recalculated Mass of Neutral Species ({final_formula_str}):")
        print(f"  C{c}: {c:2d} * {MASS_C:.5f} = {mass_C_part:9.5f} u")
        print(f"  H{h}: {h:2d} * {MASS_H:.5f} = {mass_H_part:9.5f} u")
        print(f"  Br{br}: {br:2d} * {MASS_79Br:.5f} = {mass_br_part:9.5f} u")
        print(f"  N{n}: {n:2d} * {MASS_N:.5f} = {mass_N_part:9.5f} u")
        print(f"  O{o}: {o:2d} * {MASS_O:.5f} = {mass_O_part:9.5f} u")
        print(f"  ---------------------------------------")
        print(f"  Total Neutral Mass   = {recalculated_mass:9.5f} u\n")
        
        print(f"Recalculated m/z of [M+H]+ = {recalculated_mass:.5f} + {MASS_PROTON:.5f} = {recalculated_mz:.5f} u")
        print(f"Observed m/z of [M+H]+      = {mz_observed:28.5f} u")
        print(f"Mass Error                  = {error_ppm:28.2f} ppm")
        
        print(f"\nThe molecular formula of the neutral species is {final_formula_str}.")
        return final_formula_str
    else:
        print("Could not find a matching formula within the specified parameters.")
        return None

# Run the function and get the final answer
final_answer = find_molecular_formula()
print(f"<<<{final_answer}>>>")
