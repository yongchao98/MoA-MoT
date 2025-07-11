def calculate_peptide_fragments():
    """
    Calculates the m/z of specific peptide fragments to identify
    those indicating lactylation on lysine.
    """
    # Monoisotopic masses of amino acid residues
    masses = {
        'A': 71.03711, 'R': 156.10111, 'D': 115.02694,
        'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'T': 101.04768, 'V': 99.06841
    }
    # Masses of chemical entities
    h2o_mass = 18.01056
    proton_mass = 1.00783
    lactyl_mass = 72.021

    peptide = "AVDLTKLIR"

    # --- Calculate m/z for lactylated y4 ion (y4*) ---
    # Fragment sequence: K*LIR
    y4_star_residues = masses['K'] + lactyl_mass + masses['L'] + masses['I'] + masses['R']
    y4_star_ion_mass = y4_star_residues + h2o_mass
    y4_star_mz = y4_star_ion_mass + proton_mass

    print("Calculation for the modified y4* ion (K*LIR):")
    print(f"Residue Mass(K*) = {masses['K']:.3f} (K) + {lactyl_mass:.3f} (Lactyl) = {masses['K'] + lactyl_mass:.3f}")
    print(f"Total Residue Mass (K*LIR) = {masses['K'] + lactyl_mass:.3f} + {masses['L']:.3f} (L) + {masses['I']:.3f} (I) + {masses['R']:.3f} (R) = {y4_star_residues:.3f} Da")
    print(f"Ion Mass = {y4_star_residues:.3f} + {h2o_mass:.3f} (H2O) = {y4_star_ion_mass:.3f} Da")
    print(f"m/z = {y4_star_ion_mass:.3f} + {proton_mass:.3f} (H+) = {y4_star_mz:.3f}")
    print("This calculated m/z of 601.404 matches the recorded value 601.392.\n")

    # --- Calculate m/z for b4 ion + H2O ---
    # This acts as contextual evidence, confirming the unmodified part of the sequence.
    # Fragment sequence: AVDL
    b4_residues = masses['A'] + masses['V'] + masses['D'] + masses['L']
    b4_ion_mz = b4_residues + proton_mass
    b4_plus_h2o_mz = b4_ion_mz + h2o_mass

    print("Calculation for the contextual b4+H2O ion (AVDL+H2O):")
    print(f"Total Residue Mass (AVDL) = {masses['A']:.3f} (A) + {masses['V']:.3f} (V) + {masses['D']:.3f} (D) + {masses['L']:.3f} (L) = {b4_residues:.3f} Da")
    print(f"b4 ion m/z = {b4_residues:.3f} + {proton_mass:.3f} (H+) = {b4_ion_mz:.3f}")
    print(f"m/z of b4+H2O = {b4_ion_mz:.3f} + {h2o_mass:.3f} (H2O) = {b4_plus_h2o_mz:.3f}")
    print("This calculated m/z of 417.235 matches the recorded value 417.223.\n")

    print("Conclusion: The pair of m/z values 601.392 (the modified y4* fragment) and 417.223 (a contextual unmodified fragment) together indicate that the lysine is lactylated.")

calculate_peptide_fragments()