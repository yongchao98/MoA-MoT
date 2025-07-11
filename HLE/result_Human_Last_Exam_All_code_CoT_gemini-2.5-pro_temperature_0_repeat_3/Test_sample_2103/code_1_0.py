import math

def find_molecular_formula():
    """
    Determines the molecular formula of a neutral species from its
    high-resolution mass spectrometry data.
    """
    # --- Step 1: Define constants and input data ---
    # Exact masses of the lightest stable isotopes
    MASS_H = 1.007825032
    MASS_C = 12.000000000
    MASS_N = 14.003074004
    MASS_O = 15.994914620
    MASS_79Br = 78.9183371

    # Given experimental data
    mz_protonated_mono = 1108.70902
    num_br = 6

    print("--- Analysis Step-by-Step ---")
    print("1. Isotopic pattern 1:6:15:20:15:6:1 indicates 6 Bromine atoms.")
    print(f"   Number of Bromine atoms (Br) = {num_br}\n")

    # --- Step 2: Calculate the mass of the neutral molecule ---
    mass_neutral_mono = mz_protonated_mono - MASS_H
    print(f"2. Calculate the mass of the neutral monoisotopic molecule [M]:")
    print(f"   Mass [M+H]+ = {mz_protonated_mono}")
    print(f"   Mass H      = {MASS_H}")
    print(f"   Mass [M]    = {mz_protonated_mono} - {MASS_H} = {mass_neutral_mono:.6f}\n")

    # --- Step 3: Determine the mass of the non-bromine portion ---
    mass_br_part = num_br * MASS_79Br
    target_mass_residual = mass_neutral_mono - mass_br_part
    print(f"3. Subtract the mass of the {num_br} Bromine atoms (⁷⁹Br):")
    print(f"   Mass of {num_br} x ⁷⁹Br = {num_br} * {MASS_79Br} = {mass_br_part:.6f}")
    print(f"   Residual Mass (C,H,N,O...) = {mass_neutral_mono:.6f} - {mass_br_part:.6f} = {target_mass_residual:.6f}\n")

    # --- Step 4: Brute-force search for the residual formula ---
    print("4. Searching for a C, H, N, O formula that matches the residual mass...")
    # Search constraints (based on typical natural products of this size)
    c_range = range(30, 45)
    n_range = range(0, 12, 2)  # Even numbers for even mass (Nitrogen Rule)
    o_range = range(0, 15)
    mass_tolerance_da = 0.003  # Corresponds to ~5 ppm for this mass

    found_formula = None

    for c in c_range:
        for n in n_range:
            for o in o_range:
                # Calculate mass of C, N, O part
                mass_c_n_o = c * MASS_C + n * MASS_N + o * MASS_O

                # Calculate required mass for H atoms
                rem_mass_for_h = target_mass_residual - mass_c_n_o

                # Estimate number of H atoms
                if rem_mass_for_h > 0:
                    h_float = rem_mass_for_h / MASS_H
                    h = int(round(h_float))

                    # Check if h is a reasonable integer
                    if abs(h_float - h) < 0.1:
                        # Candidate formula found, check its mass precisely
                        calculated_mass = mass_c_n_o + h * MASS_H
                        mass_diff = abs(calculated_mass - target_mass_residual)

                        if mass_diff < mass_tolerance_da:
                            # Check for valid Double Bond Equivalent (DBE)
                            # DBE = C - H/2 + N/2 + 1
                            dbe = c - h / 2.0 + n / 2.0 + 1
                            if dbe >= 0 and dbe == int(dbe):
                                found_formula = {'C': c, 'H': h, 'Br': num_br, 'N': n, 'O': o}
                                ppm_error = (mass_diff / target_mass_residual) * 1e6
                                print("\n   Match found!")
                                print(f"   Residual Formula: C{c}H{h}N{n}O{o}")
                                print(f"   Calculated Mass:  {calculated_mass:.6f}")
                                print(f"   Mass Error:       {mass_diff:.6f} Da ({ppm_error:.2f} ppm)")
                                print(f"   DBE:              {dbe}")
                                break
            if found_formula:
                break
        if found_formula:
            break
    
    # --- Step 5: Report the final result ---
    print("\n--- Final Result ---")
    if found_formula:
        c, h, br, n, o = found_formula.values()
        print("The molecular formula of the neutral species is:")
        # The final equation showing each number
        print(f"C{c} H{h} Br{br} N{n} O{o}")
        return f"C{c}H{h}Br{br}N{n}O{o}"
    else:
        print("No suitable molecular formula was found within the given constraints.")
        return None

# Run the function to find the answer
final_formula = find_molecular_formula()
print(f"\n<<<C{final_formula[1:]}>>>")
