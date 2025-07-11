import math

def find_molecular_formula():
    """
    Calculates the molecular formula of a neutral species based on high-resolution
    mass spectrometry data of its protonated ion.
    """
    # Step 1: Define constants and input data
    # High-precision atomic masses of the lightest stable isotopes
    MASS_C = 12.0000000
    MASS_H = 1.0078250
    MASS_N = 14.0030740
    MASS_O = 15.9949146
    MASS_BR79 = 78.9183371

    # Input from the problem description
    protonated_ion_mz = 1108.70902
    num_bromine = 6
    
    print("--- Step 1: Calculation Setup ---")
    print(f"Isotopic pattern indicates {num_bromine} Bromine atoms.")
    print(f"Mass of protonated monoisotopic ion [M+H]+: {protonated_ion_mz}")
    print(f"Mass of H: {MASS_H}")
    print(f"Mass of 79Br: {MASS_BR79}")
    print("-" * 35 + "\n")

    # Step 2: Calculate the mass of the neutral monoisotopic molecule (M)
    neutral_mono_mass = protonated_ion_mz - MASS_H
    
    print("--- Step 2: Calculate Neutral Mass (M) ---")
    print(f"{protonated_ion_mz} ([M+H]+) - {MASS_H} (H) = {neutral_mono_mass:.7f} (M)")
    print("-" * 35 + "\n")
    
    # Step 3: Calculate the residual mass (the CxHyNzOw part)
    mass_of_bromines = num_bromine * MASS_BR79
    residual_mass = neutral_mono_mass - mass_of_bromines

    print("--- Step 3: Calculate Residual Mass (CxHyNzOw part) ---")
    print(f"Mass of 6 * 79Br = 6 * {MASS_BR79} = {mass_of_bromines:.7f}")
    print(f"Residual Mass = {neutral_mono_mass:.7f} - {mass_of_bromines:.7f} = {residual_mass:.7f}")
    print("-" * 35 + "\n")

    # Step 4 & 5: Search for formulas and apply chemical rules
    print("--- Step 4 & 5: Searching for plausible formulas ---")
    
    # Search parameters
    tolerance_ppm = 3.0
    mass_tolerance_da = (tolerance_ppm / 1_000_000) * residual_mass
    
    # Define search limits for C, N, O atoms
    # These ranges are generous for a residual mass of ~634 Da
    max_c = 40
    max_n = 10
    max_o = 15
    
    possible_formulas = []

    # Nominal mass of the neutral molecule is 1107 (odd) -> N count must be odd
    for n in range(1, max_n, 2):
        for o in range(max_o):
            for c in range(1, max_c):
                
                # Calculate remaining mass for H
                current_mass = c * MASS_C + n * MASS_N + o * MASS_O
                mass_for_h = residual_mass - current_mass
                if mass_for_h < 0:
                    break

                # Calculate number of H atoms. Must be odd for RDBE to be integer
                h = round(mass_for_h / MASS_H)
                if h <= 0 or h % 2 == 0:
                    continue

                # Check if mass is within tolerance
                candidate_mass = current_mass + h * MASS_H
                if abs(candidate_mass - residual_mass) <= mass_tolerance_da:
                    
                    # Calculate Rings plus Double Bonds Equivalency (RDBE)
                    # RDBE = C - H/2 - X/2 + N/2 + 1
                    rdb = c - h/2 - num_bromine/2 + n/2 + 1
                    
                    # RDBE must be a non-negative integer
                    if rdb >= 0 and rdb.is_integer():
                        formula = f"C{c}H{h}N{n}O{o}Br{num_bromine}"
                        error = candidate_mass - residual_mass
                        error_ppm = (error / residual_mass) * 1e6
                        
                        possible_formulas.append({
                            "formula": formula,
                            "rdb": int(rdb),
                            "error_ppm": error_ppm
                        })

    if not possible_formulas:
        print("No plausible formula found with the given constraints.")
    else:
        print("Found the following plausible formula(s):")
        for p in possible_formulas:
            print(f"  > Formula: {p['formula']}")
            print(f"    RDBE value: {p['rdb']}")
            print(f"    Mass Error: {p['error_ppm']:.3f} ppm")
            
            # Since only one result is expected, we print it out as the final equation
            print("\n--- Final Equation ---")
            print(f"Molecular Formula of Neutral Species: {p['formula']}")
            print(f"Final Calculated Mass: {candidate_mass + mass_of_bromines:.7f}")
            print(f"Difference from target: {(candidate_mass + mass_of_bromines) - neutral_mono_mass:.7f} Da")


find_molecular_formula()