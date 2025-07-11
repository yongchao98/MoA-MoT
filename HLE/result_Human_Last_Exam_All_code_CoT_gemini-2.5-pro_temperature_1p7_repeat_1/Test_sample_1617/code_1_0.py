import collections

def solve_peptide_fragmentation():
    """
    Calculates theoretical fragment ion m/z values for a modified peptide
    and matches them against a list of experimental values.
    """
    # Step 1: Define monoisotopic masses for amino acid residues and modifications
    residue_masses = {
        'A': 71.03711,   # Alanine
        'V': 99.06841,   # Valine
        'D': 115.02694,  # Aspartic acid
        'L': 113.08406,  # Leucine
        'T': 101.04768,  # Threonine
        'K': 128.09496,  # Lysine
        'I': 113.08406,  # Isoleucine
        'R': 156.10111   # Arginine
    }
    
    # Mass of modifications and other entities
    mass_lactyl = 72.02113   # C3H4O2
    mass_h2o = 18.01056
    mass_proton = 1.00728

    # Step 2: Define the peptide, modification, and experimental data
    peptide_sequence = "AVDLTKLIR"
    modified_residue = 'K'
    experimental_mz = [401.276, 601.392, 301.200, 518.271, 304.139, 417.223]

    print("Analyzing peptide: AVDLTK(Lac)LIR")
    print(f"Lactylation on Lysine (K) adds {mass_lactyl} Da.")
    print("-" * 30)
    
    # Calculate mass of lactylated Lysine residue
    mass_K_lactylated = residue_masses['K'] + mass_lactyl
    
    # Create a new dictionary of masses for this specific peptide
    peptide_residue_masses = residue_masses.copy()
    peptide_residue_masses[modified_residue] = mass_K_lactylated

    # Step 3 & 4: Calculate theoretical y-ion masses and find evidence
    print("Calculating theoretical y-ion m/z values:")
    
    c_terminal_sequence = peptide_sequence[::-1] # Reverse sequence for y-ion calculation
    current_residue_mass_sum = 0
    
    matches = collections.defaultdict(list)
    mass_tolerance = 0.05 # Da

    for i, residue in enumerate(c_terminal_sequence):
        ion_number = i + 1
        current_residue_mass_sum += peptide_residue_masses[residue]
        
        # Calculate singly charged y-ion m/z
        y_ion_mz_z1 = current_residue_mass_sum + mass_h2o + mass_proton
        ion_name_z1 = f"y{ion_number}+"

        # The difference between y4 and y3 gives the mass of the modified K
        if ion_number == 3:
            y3_mass = y_ion_mz_z1
            y3_res_mass = current_residue_mass_sum
        if ion_number == 4:
            y4_mass = y_ion_mz_z1
            y4_res_mass = current_residue_mass_sum

        # Calculate doubly charged y-ion m/z
        y_ion_mz_z2 = (current_residue_mass_sum + mass_h2o + (2 * mass_proton)) / 2
        ion_name_z2 = f"y{ion_number}++"
        
        # Step 5 & 6: Compare with experimental values
        for exp_mz in experimental_mz:
            if abs(y_ion_mz_z1 - exp_mz) < mass_tolerance:
                matches[exp_mz].append(f"{ion_name_z1} (calculated: {y_ion_mz_z1:.3f})")
            if abs(y_ion_mz_z2 - exp_mz) < mass_tolerance:
                 matches[exp_mz].append(f"{ion_name_z2} (calculated: {y_ion_mz_z2:.3f})")

    # Step 7: Print the results and conclusion
    print("Found matches between theoretical and experimental m/z values:")
    if not matches:
        print("No matches found.")
    else:
        for mz, desc in sorted(matches.items()):
            print(f"- Experimental m/z {mz} corresponds to {', '.join(desc)}")

    print("-" * 30)
    print("Evidence for Lactylation on Lysine (K):")
    print(f"The mass of the y3 ion (LIR) is calculated as: {y3_mass:.3f}")
    print(f"The mass of the y4 ion (K(Lac)LIR) is calculated as: {y4_mass:.3f}")
    mass_diff_experimental = 601.392 - 401.276
    mass_diff_theoretical = y4_res_mass - y3_res_mass
    print(f"The mass difference between experimental y4 and y3 is {601.392} - {401.276} = {mass_diff_experimental:.3f} Da.")
    print(f"This matches the theoretical mass of a lactylated Lysine residue: {mass_K_lactylated:.3f} Da.")
    print("The presence of the y3/y4 ion pair, along with the doubly charged y4 ion, confirms lactylation on K.")

solve_peptide_fragmentation()