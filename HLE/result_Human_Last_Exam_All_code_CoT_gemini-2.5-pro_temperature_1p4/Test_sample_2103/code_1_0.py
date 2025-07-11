import math

def find_molecular_formula():
    """
    Calculates a molecular formula based on high-resolution mass spectrometry data.
    """
    # --- Step 0: Define Constants and Input Data ---
    # Monoisotopic masses of elements (in Dalton)
    C_MASS = 12.000000
    H_MASS = 1.007825
    N_MASS = 14.003074
    O_MASS = 15.994915
    BR_MASS = 78.918337

    # Input data from the problem
    mz_observed = 1108.70902  # m/z of the monoisotopic [M+H]+ ion
    num_br = 6 # Determined from the isotopic pattern 1:6:15:20:15:6:1

    print("--- Step-by-step Analysis ---")

    # --- Step 1: Calculate the mass of the neutral molecule ---
    mass_neutral_total = mz_observed - H_MASS
    print(f"1. Observed m/z of [M+H]+: {mz_observed}")
    print(f"   Subtracting mass of H ({H_MASS:.6f}):")
    print(f"   Calculated monoisotopic mass of the neutral molecule M: {mass_neutral_total:.6f} Da")
    print("-" * 30)

    # --- Step 2: Calculate the mass of the non-bromine part of the molecule ---
    mass_of_br6 = num_br * BR_MASS
    mass_residual_target = mass_neutral_total - mass_of_br6
    print(f"2. Number of bromine atoms from isotopic pattern: {num_br}")
    print(f"   Subtracting mass of {num_br} x Br ({mass_of_br6:.6f}):")
    print(f"   Target residual mass for the CxHyNzO... part: {mass_residual_target:.6f} Da")
    print("-" * 30)

    # --- Step 3: Search for the CxHyNzO.. formula ---
    print("3. Searching for a matching C, H, N, O formula...")
    # Set a realistic mass tolerance for high-resolution MS (e.g., 5 ppm)
    tolerance_ppm = 5.0
    tolerance_da = (tolerance_ppm / 1_000_000) * mass_residual_target

    # Set reasonable search ranges for the number of atoms
    # Ranges are chosen based on the residual mass and chemical intuition for natural products
    max_c = int(mass_residual_target / C_MASS) + 1
    max_n = 20
    max_o = 20
    
    solution_found = False
    final_formula = {}

    # Iterate from heaviest elements to prune search space faster
    for o_count in range(max_o + 1):
        mass_o = o_count * O_MASS
        for n_count in range(max_n + 1):
            mass_on = mass_o + n_count * N_MASS
            for c_count in range(max_c + 1):
                mass_onc = mass_on + c_count * C_MASS
                
                # Calculate required mass for H and the number of H atoms
                mass_for_h = mass_residual_target - mass_onc
                if mass_for_h < 0:
                    continue # No space left for hydrogen
                
                h_count = round(mass_for_h / H_MASS)
                
                # Calculate the mass of the candidate formula part
                candidate_mass = mass_onc + h_count * H_MASS
                
                # Check if the mass is within the tolerance
                if abs(candidate_mass - mass_residual_target) <= tolerance_da:
                    # Check chemical plausibility using Degree of Unsaturation (DBE)
                    # DBE = C - (H+X)/2 + N/2 + 1, where X is the number of halogens
                    dbe = c_count - (h_count + num_br) / 2.0 + n_count / 2.0 + 1
                    
                    # DBE must be a non-negative integer for a stable neutral molecule
                    if dbe >= 0 and abs(dbe - round(dbe)) < 0.001:
                        final_formula = {'C': c_count, 'H': h_count, 'N': n_count, 'O': o_count, 'Br': num_br}
                        solution_found = True
                        break
            if solution_found:
                break
        if solution_found:
            break

    # --- Step 4: Display the final result ---
    print("   Search complete.")
    print("-" * 30)
    if solution_found:
        c, h, n, o, br = final_formula['C'], final_formula['H'], final_formula['N'], final_formula['O'], final_formula['Br']
        print(f"4. Found a plausible formula: C{c}H{h}N{n}O{o}")
        print("\n--- FINAL RESULT ---")
        full_formula = f"C{c}H{h}N{n}O{o}Br{br}"
        print(f"The molecular formula of the neutral species is: {full_formula}")

        # Final equation as requested
        print("\nVerification of mass calculation:")
        total_calc_mass = (c*C_MASS) + (h*H_MASS) + (n*N_MASS) + (o*O_MASS) + (br*BR_MASS)
        print(f"({c} * {C_MASS}) + ({h} * {H_MASS:.6f}) + ({n} * {N_MASS:.6f}) + ({o} * {O_MASS:.6f}) + ({br} * {BR_MASS:.6f}) = {total_calc_mass:.6f} Da")
        print(f"This matches the experimentally derived neutral mass ({mass_neutral_total:.6f} Da).")
        print(f"\n<<<C{c}H{h}N{n}O{o}Br{br}>>>")
    else:
        print("No plausible formula was found within the given parameters.")

if __name__ == '__main__':
    find_molecular_formula()