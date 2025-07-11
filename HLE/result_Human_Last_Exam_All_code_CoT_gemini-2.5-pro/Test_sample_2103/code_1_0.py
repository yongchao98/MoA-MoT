import math

def find_molecular_formula():
    """
    Calculates the molecular formula of a neutral species based on high-resolution
    mass spectrometry data of its protonated form.
    """
    # --- Provided Data ---
    m_z_observed = 1108.70902
    
    # --- Step 1: Define Isotopic Masses and Constants ---
    # The isotopic pattern 1:6:15:20:15:6:1 at +2 Da increments indicates 6 Bromine atoms.
    num_Br = 6
    
    m_H = 1.00782503223      # Mass of Hydrogen atom
    m_C = 12.0000000         # Mass of Carbon-12
    m_N = 14.00307400443     # Mass of Nitrogen-14
    m_O = 15.994914619       # Mass of Oxygen-16
    m_Br79 = 78.9183371      # Mass of Bromine-79

    # --- Step 2: Calculate Mass of the Neutral Molecule ---
    mass_neutral_mono = m_z_observed - m_H
    print(f"The task is to find the molecular formula for a neutral molecule.")
    print(f"The observed m/z of the protonated ion [M+H]+ is {m_z_observed:.5f} Da.")
    print(f"The calculated monoisotopic mass of the neutral species M is {mass_neutral_mono:.5f} Da.\n")

    # --- Step 3: Determine the Mass of the Unknown Portion ---
    mass_of_br6 = num_Br * m_Br79
    mass_remainder = mass_neutral_mono - mass_of_br6
    print(f"Based on the isotopic pattern, the molecule contains {num_Br} Bromine atoms.")
    print(f"The mass of the non-Bromine portion (CxHyNzOw) is {mass_remainder:.5f} Da.\n")
    print("Searching for a valid chemical formula that matches this mass...")

    # --- Step 4 & 5: Computational Search with Chemical Rules ---
    found = False
    tolerance_da = 0.003  # 3 mDa tolerance, common for HRMS

    # Define reasonable search ranges for a natural product of this size
    c_min, c_max = 25, 50
    n_min, n_max = 1, 15  # Odd numbers only, per Nitrogen Rule for an odd total mass
    o_min, o_max = 5, 20

    for num_n in range(n_min, n_max, 2):
        for num_o in range(o_min, o_max):
            for num_c in range(c_min, c_max):
                # Calculate mass accounted for by C, N, O
                mass_c_n_o = num_c * m_C + num_n * m_N + num_o * m_O
                mass_for_h = mass_remainder - mass_c_n_o
                
                if mass_for_h < 0:
                    continue
                
                # Calculate required number of H atoms
                num_h = int(round(mass_for_h / m_H))

                # Validate H count
                if num_h <= 0:
                    continue

                # Parity Rule for integer DoU: If N is odd, H must be odd.
                if num_h % 2 == 0:
                    continue

                # --- Step 6: Verify the Candidate Formula ---
                calculated_mass_remainder = mass_c_n_o + num_h * m_H
                
                if abs(calculated_mass_remainder - mass_remainder) < tolerance_da:
                    # Degree of Unsaturation must be a non-negative integer
                    dou = num_c - num_h / 2 + num_n / 2 + 1
                    
                    if dou >= 0 and dou == int(dou):
                        formula_str = f"C{num_c}H{num_h}N{num_n}O{num_o}Br{num_Br}"
                        total_mass = calculated_mass_remainder + mass_of_br6
                        mass_error_ppm = ((total_mass - mass_neutral_mono) / mass_neutral_mono) * 1e6
                        
                        print(f"\nFound a plausible molecular formula: {formula_str}")
                        print(f"Let's verify the mass calculation:")
                        print(f"  Mass of {num_c} Carbon atoms:     {num_c * m_C:.5f}")
                        print(f"  Mass of {num_h} Hydrogen atoms:   {num_h * m_H:.5f}")
                        print(f"  Mass of {num_n} Nitrogen atoms:   {num_n * m_N:.5f}")
                        print(f"  Mass of {num_o} Oxygen atoms:     {num_o * m_O:.5f}")
                        print(f"  Mass of {num_Br} Bromine atoms:    {mass_of_br6:.5f}")
                        print(f"--------------------------------------------------")
                        print(f"Total calculated mass of neutral species: {total_mass:.5f} Da")
                        print(f"Target mass of neutral species:           {mass_neutral_mono:.5f} Da")
                        print(f"Mass Error: {mass_error_ppm:.2f} ppm")
                        print(f"Degree of Unsaturation (DoU): {int(dou)}")
                        
                        final_formula = f"C{num_c}H{num_h}N{num_n}O{num_o}Br{num_Br}"
                        print(f"\nThe final molecular formula is presented below.")
                        print(f"<<<{final_formula}>>>")
                        found = True
                        return

    if not found:
        print("\nNo plausible formula found with the given constraints.")

if __name__ == "__main__":
    find_molecular_formula()