import math

def find_molecular_formula():
    """
    Calculates the molecular formula of a neutral species based on high-resolution
    mass spectrometry data of its protonated ion.
    """
    # --- Part 1: Define Constants and Knowns ---

    # High precision atomic masses (monoisotopic)
    MASS_H = 1.007825
    MASS_C = 12.000000
    MASS_N = 14.003074
    MASS_O = 15.994915
    MASS_Br79 = 78.918337

    # Experimental data from the problem
    MASS_PROTONATED = 1108.70902
    NUM_Br = 6  # Deduced from the 1:6:15:20:15:6:1 isotopic pattern

    # Set a realistic mass tolerance for high-resolution MS (5 ppm)
    mass_tolerance_ppm = 5.0

    # --- Part 2: Calculate Target Masses ---

    # The mass of the neutral molecule is the mass of the protonated ion minus the mass of a proton.
    mass_neutral_target = MASS_PROTONATED - MASS_H

    # The mass tolerance in absolute units (Daltons)
    mass_tolerance_da = (mass_neutral_target * mass_tolerance_ppm) / 1_000_000

    # Subtract the mass of the 6 bromine atoms to find the mass of the C,H,N,O portion.
    mass_of_bromines = NUM_Br * MASS_Br79
    mass_remainder_target = mass_neutral_target - mass_of_bromines

    # --- Part 3: Systematically Search for the Formula ---

    solution_found = False
    final_formula = {}

    # Iterate through plausible ranges for C, N, and O atoms.
    # The nominal mass of the neutral molecule (1107) is odd, so N must be odd.
    for num_c in range(30, 40):
        for num_n in range(1, 10, 2):  # Step by 2 to only check odd numbers
            for num_o in range(4, 12):
                # Calculate the mass of the current C,N,O combination
                mass_cno = (num_c * MASS_C) + (num_n * MASS_N) + (num_o * MASS_O)
                
                # Calculate the remaining mass that must be from Hydrogen atoms
                mass_h_remaining = mass_remainder_target - mass_cno
                
                if mass_h_remaining > 0:
                    # Calculate the potential number of hydrogens (as a float)
                    num_h_float = mass_h_remaining / MASS_H
                    
                    # Check if the number of hydrogens is very close to a whole number
                    if abs(num_h_float - round(num_h_float)) < 0.01:
                        num_h = int(round(num_h_float))
                        
                        # Recalculate the mass with the integer H count to check against tolerance
                        calculated_remainder_mass = mass_cno + (num_h * MASS_H)
                        
                        if abs(calculated_remainder_mass - mass_remainder_target) < mass_tolerance_da:
                            solution_found = True
                            final_formula = {'C': num_c, 'H': num_h, 'N': num_n, 'O': num_o, 'Br': NUM_Br}
                            break
            if solution_found:
                break
        if solution_found:
            break

    # --- Part 4: Print the Final Result ---

    if solution_found:
        c, h, n, o, br = final_formula['C'], final_formula['H'], final_formula['N'], final_formula['O'], final_formula['Br']
        
        # Calculate the masses of each component for the final report
        mass_calc_C = c * MASS_C
        mass_calc_H = h * MASS_H
        mass_calc_N = n * MASS_N
        mass_calc_O = o * MASS_O
        mass_calc_Br = br * MASS_Br79
        total_calculated_mass = mass_calc_C + mass_calc_H + mass_calc_N + mass_calc_O + mass_calc_Br

        formula_string = f"C{c}H{h}N{n}O{o}Br{br}"
        error_da = total_calculated_mass - mass_neutral_target
        error_ppm = (error_da / mass_neutral_target) * 1_000_000

        print(f"The observed m/z of [M+H]⁺ is {MASS_PROTONATED:.5f}.")
        print(f"The calculated mass of the neutral species (M) is {mass_neutral_target:.5f} Da.\n")
        print(f"A plausible molecular formula that fits this mass is {formula_string}.")
        print("The mass breakdown is as follows:\n")
        
        print(f"{c:>4} C  = {c:>2} * {MASS_C:<10.6f} = {mass_calc_C:>9.5f}")
        print(f"{h:>4} H  = {h:>2} * {MASS_H:<10.6f} = {mass_calc_H:>9.5f}")
        print(f"{n:>4} N  = {n:>2} * {MASS_N:<10.6f} = {mass_calc_N:>9.5f}")
        print(f"{o:>4} O  = {o:>2} * {MASS_O:<10.6f} = {mass_calc_O:>9.5f}")
        print(f"{br:>4} Br = {br:>2} * {MASS_Br79:<10.6f} = {mass_calc_Br:>9.5f}")
        print("-" * 40)
        print(f"Total Calculated Mass for {formula_string} = {total_calculated_mass:>9.5f}\n")
        print(f"The mass error is {error_da:.5f} Da ({error_ppm:.2f} ppm), which is well within acceptable limits.")

    else:
        print("No molecular formula was found within the specified search parameters and tolerance.")

if __name__ == '__main__':
    find_molecular_formula()