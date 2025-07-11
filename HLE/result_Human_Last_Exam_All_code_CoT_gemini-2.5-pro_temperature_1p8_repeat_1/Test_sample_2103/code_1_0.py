import math

def find_molecular_formula():
    """
    Analyzes mass spectrometry data to determine the molecular formula of a neutral species.
    """
    # Step 1: Define constants and experimental data
    m_H = 1.0078250322      # Mass of Hydrogen
    m_C = 12.0000000000      # Mass of Carbon-12
    m_N = 14.0030740044      # Mass of Nitrogen-14
    m_O = 15.9949146196      # Mass of Oxygen-16
    m_Br79 = 78.9183371        # Mass of Bromine-79
    m_H_plus = 1.0072764669    # Mass of a proton

    # Isotopic signature suggests 6 bromine atoms
    num_Br = 6

    # Observed m/z for the monoisotopic, protonated species [M+H]+
    mz_observed = 1108.70902

    # Step 2: Calculate the exact mass of the neutral molecule
    mass_neutral_mono = mz_observed - m_H_plus

    # Step 3: Subtract the mass of the bromine atoms to find the remainder mass
    mass_br_total = num_Br * m_Br79
    mass_remainder_target = mass_neutral_mono - mass_br_total

    # Step 4: Search for a molecular formula (CxHyNzOa) that matches the remainder mass
    # Rules for search: odd number of Nitrogens, odd number of Hydrogens.
    ppm_tolerance = 5.0  # Tolerance for high-resolution MS
    mass_tolerance = mass_neutral_mono * ppm_tolerance / 1e6

    found_formula = None
    final_calculated_mass = 0

    # Search range for atoms, typical for this class of natural products
    for C in range(30, 40):
        for N in range(1, 10, 2):  # Odd N for odd nominal mass
            for O in range(5, 15):
                # Back-calculate the required mass for Hydrogen atoms
                mass_CNO = C * m_C + N * m_N + O * m_O
                mass_for_H = mass_remainder_target - mass_CNO

                if mass_for_H > 0:
                    num_H_float = mass_for_H / m_H
                    num_H = round(num_H_float)

                    # For DBE to be integer, if N is odd, H must also be odd
                    if num_H % 2 != 1:
                        continue

                    calculated_mass_neutral = (C * m_C + num_H * m_H + N * m_N + O * m_O) + mass_br_total

                    if abs(calculated_mass_neutral - mass_neutral_mono) <= mass_tolerance:
                        found_formula = {'C': C, 'H': num_H, 'Br': num_Br, 'N': N, 'O': O}
                        final_calculated_mass = calculated_mass_neutral
                        break
            if found_formula: break
        if found_formula: break

    # Step 5: Output the result
    print("--- Analysis of UPLC-MS Data ---")
    print(f"Observed [M+H]+ m/z: {mz_observed}")
    print(f"Calculated Neutral Monoisotopic Mass: {mass_neutral_mono:.6f} Da")
    print(f"Target mass for the non-bromine part (CxHyNzOa): {mass_remainder_target:.6f} Da")
    print("-" * 35)

    if found_formula:
        print("A matching molecular formula was found!")
        C, H, Br, N, O = found_formula['C'], found_formula['H'], found_formula['Br'], found_formula['N'], found_formula['O']
        print(f"\nMolecular Formula of Neutral Species: C{C}H{H}Br{Br}N{N}O{O}")

        # Verification part
        print("\nVerification of mass calculation:")
        m_C_tot = C * m_C
        m_H_tot = H * m_H
        m_N_tot = N * m_N
        m_O_tot = O * m_O
        print(f"Mass of {C} Carbons:    {m_C_tot: >12.6f} Da")
        print(f"Mass of {H} Hydrogens:   {m_H_tot: >12.6f} Da")
        print(f"Mass of {Br} Bromines:   {mass_br_total: >12.6f} Da")
        print(f"Mass of {N} Nitrogens:   {m_N_tot: >12.6f} Da")
        print(f"Mass of {O} Oxygens:     {m_O_tot: >12.6f} Da")
        print("------------------------------------------ +")
        print(f"Calculated Total Mass: {final_calculated_mass: >12.6f} Da")

        error_ppm = ((final_calculated_mass - mass_neutral_mono) / mass_neutral_mono) * 1e6
        print(f"\nThis calculated mass is within {error_ppm:.2f} ppm of the experimental mass.")
    else:
        print("No matching molecular formula was found within the specified search parameters and tolerance.")

if __name__ == '__main__':
    find_molecular_formula()