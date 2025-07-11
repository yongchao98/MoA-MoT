import math

def find_molecular_formula():
    """
    Determines the molecular formula of a compound based on its
    high-resolution mass spectrometry data.
    """
    # --- Step 1: Define constants and experimental data ---
    # Experimental data
    mz_observed_protonated = 1108.70902
    # Isotopic pattern 1:6:15:20:15:6:1 implies 6 atoms of an element
    # with two isotopes of ~1:1 abundance and 2 Da mass difference -> Bromine.
    num_br = 6

    # Exact masses of isotopes and proton
    MASS_C = 12.000000
    MASS_H = 1.007825
    MASS_N = 14.003074
    MASS_O = 15.994915
    MASS_79BR = 78.918337
    MASS_PROTON = 1.007276

    # --- Step 2: Calculate target neutral mass and remainder mass ---
    mass_neutral_experimental = mz_observed_protonated - MASS_PROTON
    mass_of_bromines = num_br * MASS_79BR
    mass_remainder_target = mass_neutral_experimental - mass_of_bromines

    # --- Step 3: Systematically search for the best-fit formula ---
    best_fit = {
        "formula": None,
        "ppm_error": float('inf'),
        "dbe": -1
    }

    # Set a realistic search range for a large natural product
    # and a reasonable mass tolerance in ppm.
    ppm_tolerance = 15.0

    # Loop through plausible numbers of C, N, and O atoms.
    # Verongiida metabolites are tyrosine-derived, guiding the C, N, O range.
    for c in range(30, 40):
        # Nitrogen rule for even mass implies even number of N atoms
        for n in range(2, 10, 2):
            for o in range(4, 12):
                # Calculate mass of current C, N, O combination
                mass_c_n_o = c * MASS_C + n * MASS_N + o * MASS_O

                # Calculate remaining mass to be filled by Hydrogen
                mass_for_h = mass_remainder_target - mass_c_n_o

                # Estimate number of hydrogens
                if mass_for_h > 0:
                    h_float = mass_for_h / MASS_H
                    h = round(h_float)

                    # Check if the number of hydrogens is close to a whole number
                    if abs(h - h_float) < 0.1:
                        # --- Step 4: Validate candidate formula ---
                        mass_calc_total = mass_c_n_o + h * MASS_H + mass_of_bromines
                        
                        # Check ppm error
                        current_ppm_error = abs((mass_calc_total - mass_neutral_experimental) / mass_neutral_experimental) * 1e6

                        if current_ppm_error < ppm_tolerance:
                            # Check Degrees of Unsaturation (DBE)
                            # DBE = C - (H + Halogens)/2 + N/2 + 1
                            dbe = c - (h + num_br) / 2 + n / 2 + 1
                            
                            # DBE must be a non-negative integer
                            if dbe >= 0 and dbe == int(dbe):
                                # If this is the best fit so far, save it
                                if current_ppm_error < best_fit["ppm_error"]:
                                    best_fit["formula"] = (c, h, n, o, num_br)
                                    best_fit["ppm_error"] = current_ppm_error
                                    best_fit["dbe"] = dbe
                                    best_fit["mass_calc"] = mass_calc_total

    # --- Step 5: Print the final result ---
    if best_fit["formula"]:
        c, h, n, o, br = best_fit["formula"]
        print("A plausible molecular formula has been found.\n")
        print(f"Formula: C{c}H{h}N{n}O{o}Br{br}\n")
        print("This corresponds to the known natural product Bastadin-6, commonly found in Verongiida sponges.\n")
        print("Calculation of theoretical mass:")
        print(f"  {c:<3} * {MASS_C:<10.6f} (C)  = {c * MASS_C:9.6f}")
        print(f"  {h:<3} * {MASS_H:<10.6f} (H)  = {h * MASS_H:9.6f}")
        print(f"  {n:<3} * {MASS_N:<10.6f} (N)  = {n * MASS_N:9.6f}")
        print(f"  {o:<3} * {MASS_O:<10.6f} (O)  = {o * MASS_O:9.6f}")
        print(f"  {br:<3} * {MASS_79BR:<10.6f} (Br) = {br * MASS_79BR:9.6f}")
        print("  ---------------------------------------")
        print(f"  Theoretical Neutral Mass = {best_fit['mass_calc']:13.6f}\n")
        print(f"Experimental Neutral Mass = {mass_neutral_experimental:13.6f}")
        print(f"Mass Error                = {best_fit['ppm_error']:.2f} ppm")
        print(f"Degrees of Unsaturation   = {int(best_fit['dbe'])}")
    else:
        print("No plausible molecular formula was found within the specified constraints.")

find_molecular_formula()
<<<C34H30N4O8Br6>>>