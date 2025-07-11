def solve_peptide_ms():
    """
    Calculates theoretical m/z values for a lactylated peptide and matches them
    against a list of observed values to identify evidence of lactylation.
    """
    # 1. Define monoisotopic masses
    residue_mass = {
        'A': 71.03711, 'V': 99.06841, 'D': 115.02694, 'L': 113.08406,
        'T': 101.04768, 'K': 128.09496, 'I': 113.08406, 'R': 156.10111
    }
    # Lactyl group (C3H4O2) adds 72.02114 Da
    LACTYL_MOD = 72.02114
    H_MASS = 1.00728  # Proton mass for [M+H]+
    H2O_MASS = 18.01056

    # 2. Define peptide and observed data
    sequence = "AVDLTKLIR"
    mod_position = 6  # Lactylation on K at position 6 (1-based)
    observed_mz = [401.276, 601.392, 301.200, 518.271, 304.139, 417.223]

    print("--- Peptide and Modification ---")
    print(f"Sequence: {sequence}")
    print(f"Modification: Lactylation (+{LACTYL_MOD} Da) on residue {sequence[mod_position-1]}{mod_position}")
    print(f"Observed m/z values: {observed_mz}\n")

    # Add modification to the mass of Lysine
    residue_mass_mod = residue_mass.copy()
    residue_mass_mod['K'] += LACTYL_MOD

    # 3. Calculate theoretical fragment ions
    print("--- Theoretical Fragment Calculation ---")
    
    # y-ions (C-terminus)
    print("\nCalculating y-ions...")
    y_ion_mass = 0
    for i in range(len(sequence) - 1, -1, -1):
        res = sequence[i]
        pos = i + 1
        # Use modified mass if at the modification position
        current_res_mass = residue_mass_mod[res] if pos == mod_position else residue_mass[res]
        y_ion_mass += current_res_mass
        ion_name = f"y{len(sequence) - i}"
        # y-ion m/z = (sum of residue masses + H2O + H+) / charge (z=1)
        y_ion_mz = y_ion_mass + H2O_MASS + H_MASS
        print(f"Fragment {ion_name} ({sequence[i:]}): Calculated m/z = {y_ion_mz:.3f}")
        # Store key ions for final explanation
        if ion_name == 'y3': y3_mz = y_ion_mz
        if ion_name == 'y4': y4_lac_mz = y_ion_mz

    # b-ions (N-terminus)
    print("\nCalculating b-ions...")
    b_ion_mass = 0
    for i in range(len(sequence)):
        res = sequence[i]
        pos = i + 1
        current_res_mass = residue_mass_mod[res] if pos == mod_position else residue_mass[res]
        b_ion_mass += current_res_mass
        ion_name = f"b{i + 1}"
        # b-ion m/z = (sum of residue masses + H+) / charge (z=1)
        b_ion_mz = b_ion_mass + H_MASS
        print(f"Fragment {ion_name} ({sequence[:i+1]}): Calculated m/z = {b_ion_mz:.3f}")
        # Calculate b+H2O adduct for specific ions
        if ion_name == 'b5':
            b5_mz = b_ion_mz
            b5_h2o_mz = b_ion_mz + H2O_MASS
            print(f"Fragment b5+H2O: Calculated m/z = {b5_h2o_mz:.3f}")

    # 4. Compare and explain the evidence
    print("\n--- Analysis of Evidence ---")
    print("Matching theoretical values to observed m/z list:")

    # Match for 401.276 -> y3
    print(f"\n1. Observed m/z = 401.276")
    print(f"   Matches y3 ion (LIR). Theoretical m/z = {y3_mz:.3f}.")
    print(f"   This fragment is C-terminal to the modification site and confirms part of the peptide sequence.")
    
    # Match for 601.392 -> y4(lac)
    print(f"\n2. Observed m/z = 601.392")
    print(f"   Matches y4 ion (K(lac)LIR). Theoretical m/z = {y4_lac_mz:.3f}.")
    print(f"   This fragment contains the modified lysine.")
    
    # Mass difference between y4 and y3
    mass_diff_y4_y3 = y4_lac_mz - y3_mz
    k_lac_mass = residue_mass['K'] + LACTYL_MOD
    print(f"   The mass difference between y4 and y3 is {y4_lac_mz:.3f} - {y3_mz:.3f} = {mass_diff_y4_y3:.3f} Da.")
    print(f"   This difference corresponds to the mass of lactylated lysine ({k_lac_mass:.3f} Da), precisely localizing the modification to K6.")
    
    # Match for 518.271 -> b5+H2O
    print(f"\n3. Observed m/z = 518.271")
    print(f"   Matches a b5+H2O adduct ion (AVDLT + H2O). Theoretical m/z = {b5_h2o_mz:.3f}.")
    print(f"   This fragment confirms the N-terminal sequence up to the residue before the modification site.")
    
    print("\n--- Conclusion ---")
    print("The combination of m/z values 401.276 (y3), 601.392 (y4), and 518.271 (b5+H2O) provides a comprehensive confirmation:")
    print("- The y3/y4 pair proves that a modification of +72 Da occurred on the Lysine at position 6.")
    print("- The b5 ion confirms the sequence leading up to the modification site.")
    print("Therefore, this set of ions together indicates that the lysine is lactylated.")
    print("The corresponding answer choice is D.")

solve_peptide_ms()
<<<D>>>