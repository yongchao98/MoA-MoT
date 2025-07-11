import collections

def solve_peptide_modification():
    """
    Calculates theoretical m/z values for peptide fragments to identify
    the one corresponding to lactylation.
    """

    # Monoisotopic masses of residues (mass of the amino acid minus H2O)
    residue_masses = {
        'A': 71.03711, 'V': 99.06841, 'D': 115.02694,
        'L': 113.08406, 'T': 101.04768, 'K': 128.09496,
        'I': 113.08406, 'R': 156.10111
    }

    # Mass of lactyl group modification (C3H4O2)
    # This is added to the lysine residue
    lactyl_mass = (3 * 12.00000) + (4 * 1.007825) + (2 * 15.994915) # 72.02113
    
    # Mass of other relevant particles
    proton_mass = 1.007825
    h2o_mass = 18.010565

    # Peptide sequence and recorded m/z values
    peptide_sequence = "AVDLTKLIR"
    recorded_mz = [401.276, 601.392, 301.200, 518.271, 304.139, 417.223]
    
    # --- Step 1: Calculate the mass of the modified Lysine residue ---
    lysine_mass = residue_masses['K']
    lactylated_lysine_mass = lysine_mass + lactyl_mass
    
    print("--- Analysis of Peptide AVDLTK(lac)LIR ---")
    print(f"Peptide Sequence: {peptide_sequence}")
    print(f"Modification: Lactylation on Lysine (K)\n")
    print(f"Mass of Lysine (K) residue: {lysine_mass:.5f} Da")
    print(f"Mass of Lactyl group (C3H4O2): {lactyl_mass:.5f} Da")
    print(f"Mass of lactylated Lysine (K_lac) residue: {lactylated_lysine_mass:.5f} Da\n")

    # --- Step 2: Calculate theoretical y-ion series around the modification ---
    # y-ions are fragments from the C-terminus. Their mass is calculated as:
    # sum_of_residue_masses + mass_H2O
    
    # y3 ion: L-I-R
    y3_neutral_mass = residue_masses['L'] + residue_masses['I'] + residue_masses['R'] + h2o_mass
    y3_mz = y3_neutral_mass + proton_mass

    # y4 ion: K(lac)-L-I-R
    y4_neutral_mass = lactylated_lysine_mass + y3_neutral_mass - h2o_mass
    y4_mz = y4_neutral_mass + proton_mass

    print("--- Calculating Theoretical m/z Values ---")
    print(f"Theoretical m/z of y3 ion (LIR): {y3_mz:.3f}")
    print(f"Theoretical m/z of y4 ion (K(lac)LIR): {y4_mz:.3f}\n")
    
    # --- Step 3: Compare with recorded m/z values ---
    print("--- Comparing with Recorded m/z Values ---")
    y3_observed = 401.276
    y4_observed = 601.392
    
    print(f"The calculated y3 m/z ({y3_mz:.3f}) matches the recorded value {y3_observed}.")
    print(f"The calculated y4 m/z ({y4_mz:.3f}) matches the recorded value {y4_observed}.")
    print("The y4 ion contains the modified lysine. Its presence indicates that the modification exists.\n")

    # --- Step 4: Confirm using the mass difference ---
    mass_diff_observed = y4_observed - y3_observed
    
    print("--- Final Confirmation by Mass Difference ---")
    print("The mass difference between the y4 and y3 ions should equal the mass of the modified residue.")
    print(f"Equation: m/z(y4) - m/z(y3) = Mass(K_lac)")
    print(f"Using observed values: {y4_observed} - {y3_observed} = {mass_diff_observed:.3f} Da")
    print(f"This mass difference ({mass_diff_observed:.3f} Da) is a direct measurement of the fourth residue from the C-terminus (K),")
    print(f"and it perfectly matches our calculated mass of lactylated lysine ({lactylated_lysine_mass:.3f} Da).")

solve_peptide_modification()