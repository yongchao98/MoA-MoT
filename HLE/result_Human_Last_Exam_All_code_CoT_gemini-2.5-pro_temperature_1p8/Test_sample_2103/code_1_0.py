import math

def find_molecular_formula():
    """
    This script determines the molecular formula of a compound based on its
    high-resolution mass spectrometry data.
    """

    # --- Define constants for lightest isotopes ---
    H_MASS = 1.007825
    C_MASS = 12.000000
    N_MASS = 14.003074
    O_MASS = 15.994915
    BR_79_MASS = 78.918337
    
    # --- Step 1: Analyze Input Data ---
    # The isotopic pattern 1:6:15:20:15:6:1 corresponds to the binomial coefficient
    # for n=6, (a+b)^6. The 2 amu spacing indicates an element with two
    # common isotopes of roughly 1:1 abundance and a mass difference of 2.
    # This signature perfectly matches Bromine (⁷⁹Br and ⁸¹Br).
    NUM_BR = 6
    
    # Observed m/z for the protonated monoisotopic peak [M+H]+
    MZ_OBSERVED = 1108.70902

    # --- Step 2: Calculate Masses ---
    # Calculate the mass of the neutral monoisotopic molecule M
    # Mass(M) = Mass([M+H]+) - Mass(H)
    mass_m_neutral = MZ_OBSERVED - H_MASS
    
    # Calculate the mass of the 6 bromine atoms (using the lightest isotope ⁷⁹Br)
    mass_br_total = NUM_BR * BR_79_MASS
    
    # Calculate the mass of the remainder of the molecule (R = C, H, N, O)
    mass_r_residual = mass_m_neutral - mass_br_total

    # --- Step 3: Apply Chemical Rules to Constrain Search ---
    # Nitrogen Rule: Applies to the neutral molecule M.
    # The nominal mass of M is floor(1107.701195) = 1107, which is ODD.
    # For a neutral molecule, an odd nominal mass implies an ODD number of nitrogen atoms.
    
    # Hydrogen Parity Rule: Derived from the nominal mass formula.
    # Nom_Mass(M) = 12*c + 1*h + 14*n + 16*o + 79*br
    # 1107 = 12*c + h + 14*n + 16*o + 6*79
    # odd = even + h + even + even + even
    # This implies that the number of hydrogen atoms 'h' must be ODD.
    
    # --- Step 4: Iterative Search for the Formula of Remainder R ---
    print("--- Molecular Formula Determination ---")
    print("Step 1: Isotopic Pattern (1:6:15:20:15:6:1) -> 6 Bromine atoms.")
    print(f"Step 2: Observed m/z [M+H]+ = {MZ_OBSERVED}")
    print(f"Step 3: Calculated neutral monoisotopic mass M = {mass_m_neutral:.6f}")
    print(f"Step 4: Calculated residual mass (M - Br6) = {mass_r_residual:.6f}")
    print("\n--- Searching for a valid formula C(c)H(h)N(n)O(o) for the residual mass ---")
    print("Constraints: N count must be odd, H count must be odd.")

    PPM_TOLERANCE = 5.0
    solutions = []

    # Set realistic loop ranges
    max_o = int(mass_r_residual / O_MASS) + 1
    max_n = int(mass_r_residual / N_MASS) + 1
    max_c = int(mass_r_residual / C_MASS) + 1
    
    # Loop with N as odd numbers only
    for num_n in range(1, max_n, 2):
        mass_from_n = num_n * N_MASS
        for num_o in range(max_o):
            mass_from_o = num_o * O_MASS
            for num_c in range(1, max_c):
                mass_from_c = num_c * C_MASS
                
                # Calculate required mass for H
                mass_for_h = mass_r_residual - (mass_from_n + mass_from_o + mass_from_c)
                if mass_for_h < 0:
                    continue

                # Calculate number of H atoms and check if it's odd
                num_h = round(mass_for_h / H_MASS)
                if num_h <= 0 or num_h % 2 != 1:
                    continue

                # Calculate theoretical mass and error
                theoretical_mass = mass_from_c + num_h * H_MASS + mass_from_n + mass_from_o
                error_ppm = abs(theoretical_mass - mass_r_residual) / mass_r_residual * 1e6
                
                if error_ppm < PPM_TOLERANCE:
                    # Final check: Degree of Unsaturation (DBE) must be a non-negative integer
                    # DBE = C - H/2 - X/2 + N/2 + 1 (X = halogens)
                    dbe = num_c - num_h / 2.0 - NUM_BR / 2.0 + num_n / 2.0 + 1
                    if dbe >= 0 and abs(dbe - round(dbe)) < 0.01:
                        solutions.append({
                            'c': num_c, 'h': num_h, 'n': num_n, 'o': num_o,
                            'dbe': dbe, 'error_ppm': error_ppm
                        })

    # --- Step 5: Display the Best Solution ---
    if not solutions:
        print("\nNo plausible formula found within the given constraints.")
        return

    best_solution = min(solutions, key=lambda x: x['error_ppm'])
    c, h, n, o = best_solution['c'], best_solution['h'], best_solution['n'], best_solution['o']
    
    print("\n--- Result ---")
    print(f"Found Formula: C{c} H{h} Br{NUM_BR} N{n} O{o}")
    print(f"Degree of Unsaturation: {best_solution['dbe']:.1f}")

    print("\n--- Mass Calculation Verification ---")
    mass_c = c * C_MASS
    mass_h = h * H_MASS
    mass_br = NUM_BR * BR_79_MASS
    mass_n = n * N_MASS
    mass_o = o * O_MASS
    
    total_neutral_mass = mass_c + mass_h + mass_br + mass_n + mass_o
    protonated_mass = total_neutral_mass + H_MASS
    error = protonated_mass - MZ_OBSERVED
    error_ppm = (error / MZ_OBSERVED) * 1e6
    
    final_equation = {
        f"C{c}": (c, C_MASS, mass_c),
        f"H{h}": (h, H_MASS, mass_h),
        f"Br{NUM_BR}": (NUM_BR, BR_79_MASS, mass_br),
        f"N{n}": (n, N_MASS, mass_n),
        f"O{o}": (o, O_MASS, mass_o),
    }

    for key, values in final_equation.items():
        print(f"{key:<5}: {values[0]:2d} * {values[1]:10.6f} = {values[2]:10.6f}")

    print("------------------------------------------")
    print(f"Calculated Neutral Mass:   {total_neutral_mass:12.6f}")
    print(f"Calculated [M+H]+ Mass:   {protonated_mass:12.6f}")
    print(f"Observed   [M+H]+ Mass:   {MZ_OBSERVED:12.6f}")
    print("==========================================")
    print(f"Mass Error:                {error:10.6f} Da ({error_ppm:.2f} ppm)")

    global final_formula_str
    final_formula_str = f"C{c}H{h}Br{NUM_BR}N{n}O{o}"


# This is to capture the final answer for the required format.
final_formula_str = "Formula not found"
find_molecular_formula()
# The required format is <<<ANSWER>>>.
# Although the print statements above give a full report, this line provides
# the final answer in the specific format requested by the prompt.
print(f"\n<<<{final_formula_str}>>>")
