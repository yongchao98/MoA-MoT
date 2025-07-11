import math

def find_molecular_formula():
    """
    Calculates and prints the molecular formula based on MS data.
    """
    # Step 0: Define constants and experimental data
    # Exact atomic masses (monoisotopic)
    M_H = 1.007825
    M_C = 12.000000
    M_N = 14.003074
    M_O = 15.994915
    M_BR79 = 78.9183371

    # Experimental data
    mz_observed = 1108.70902
    num_bromines = 6
    
    print("--- Task: Determine Molecular Formula from MS Data ---")
    print(f"\nStep 1: Analyze Isotopic Pattern")
    print("The 1:6:15:20:15:6:1 pattern with 2 amu spacing indicates 6 Bromine atoms.")
    print(f"Number of Bromine atoms (Br) = {num_bromines}")

    # Step 1: Calculate the mass of the non-halogen, non-proton fragment
    # For a protonated ion [M+H]+, the mass of the neutral M is mz - mass(H)
    # when using neutral atomic masses throughout.
    mass_neutral_monoisotopic = mz_observed - M_H
    mass_of_bromines = num_bromines * M_BR79
    target_mass_fragment = mass_neutral_monoisotopic - mass_of_bromines
    
    print(f"\nStep 2: Calculate Target Mass of the Remainder (CxHyNzOw...)")
    print(f"Observed m/z of [M+H]+ = {mz_observed}")
    print(f"Subtracting mass of a proton (H atom) = {M_H}")
    print(f"Mass of neutral molecule M (lightest isotopes) = {mass_neutral_monoisotopic:.5f}")
    print(f"Subtracting mass of {num_bromines} x 79Br atoms ({num_bromines} * {M_BR79:.5f}) = {mass_of_bromines:.5f}")
    print(f"Target mass for the CxHyNzOw... fragment = {target_mass_fragment:.5f}")
    
    # Step 2: Iterate through plausible element counts to find the formula
    tolerance_ppm = 5.0
    tolerance_amu = target_mass_fragment * tolerance_ppm / 1e6
    best_hit = None
    lowest_error = float('inf')

    # Heuristic search ranges
    # Nitrogen Rule: Nominal mass of neutral molecule is 1107 (odd), so N count should be odd.
    n_range = range(1, 10, 2)
    c_range = range(25, 45)
    o_range = range(5, 15)

    for n_c in c_range:
        for n_n in n_range:
            for n_o in o_range:
                # Calculate remaining mass for hydrogens
                cno_mass = n_c * M_C + n_n * M_N + n_o * M_O
                h_mass_needed = target_mass_fragment - cno_mass

                if h_mass_needed > 0:
                    # Number of H must be a positive integer
                    n_h = int(round(h_mass_needed / M_H))
                    if n_h <= 0:
                        continue
                        
                    # Calculate mass of the candidate fragment and its error
                    calculated_fragment_mass = cno_mass + n_h * M_H
                    mass_error_amu = abs(calculated_fragment_mass - target_mass_fragment)

                    if mass_error_amu < tolerance_amu and mass_error_amu < lowest_error:
                        # Candidate found, check Double Bond Equivalent (DBE)
                        # DBE = C - (H+X)/2 + N/2 + 1
                        DBE = n_c - (n_h + num_bromines) / 2.0 + n_n / 2.0 + 1
                        # If DBE is a non-negative integer or half-integer, we consider it.
                        if DBE >= 0:
                            lowest_error = mass_error_amu
                            best_hit = {'C': n_c, 'H': n_h, 'N': n_n, 'O': n_o, 
                                        'Br': num_bromines, 'DBE': DBE}
    
    print(f"\nStep 3: Search for the Best-Fit Formula")
    if best_hit:
        formula_str = f"C{best_hit['C']}H{best_hit['H']}Br{best_hit['Br']}N{best_hit['N']}O{best_hit['O']}"
        print(f"Best matching formula found: {formula_str}")
        print(f"This formula has a Double Bond Equivalent (DBE) of {best_hit['DBE']}.")
        if best_hit['DBE'] != int(best_hit['DBE']):
            print("Note: A half-integer DBE implies the neutral species is a radical.")

        # Step 4: Verify the result
        print("\nStep 4: Verify the Mass of the Proposed Formula")
        c, h, n, o, br = best_hit['C'], best_hit['H'], best_hit['N'], best_hit['O'], best_hit['Br']
        
        calc_mass_c = c * M_C
        calc_mass_h = h * M_H
        calc_mass_n = n * M_N
        calc_mass_o = o * M_O
        calc_mass_br = br * M_BR79
        
        calc_total_neutral_mass = calc_mass_c + calc_mass_h + calc_mass_n + calc_mass_o + calc_mass_br
        calc_mz = calc_total_neutral_mass + M_H

        final_error_ppm = (abs(calc_mz - mz_observed) / mz_observed) * 1e6

        print(f"Final proposed neutral formula: {formula_str}")
        print("\nMass Calculation:")
        print(f"  Mass C{c} = {c} * {M_C:.6f} = {calc_mass_c:.6f}")
        print(f"  Mass H{h} = {h} * {M_H:.6f} = {calc_mass_h:.6f}")
        print(f"  Mass Br{br} = {br} * {M_BR79:.6f} = {calc_mass_br:.6f}")
        print(f"  Mass N{n} = {n} * {M_N:.6f} = {calc_mass_n:.6f}")
        print(f"  Mass O{o} = {o} * {M_O:.6f} = {calc_mass_o:.6f}")
        print(f"  --------------------------------------------------")
        print(f"  Total Mass of Neutral Molecule M = {calc_total_neutral_mass:.6f}")
        print(f"\nCalculated m/z of [M+H]+:")
        print(f"  Mass(M) + Mass(H) = {calc_total_neutral_mass:.6f} + {M_H:.6f} = {calc_mz:.6f}")
        print(f"Experimental m/z of [M+H]+: {mz_observed:.6f}")
        print(f"Mass Error: {final_error_ppm:.2f} ppm")

    else:
        print("No suitable formula found within the given constraints and search space.")

# Execute the function to find and print the answer
find_molecular_formula()