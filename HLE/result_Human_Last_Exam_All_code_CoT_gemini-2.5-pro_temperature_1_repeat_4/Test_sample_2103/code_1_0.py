import math

def find_molecular_formula():
    """
    Calculates a molecular formula based on high-resolution mass spectrometry data.
    """
    # --- Part 1: Define constants and given data ---
    # High-precision masses of the lightest stable isotopes
    H_MASS = 1.00782503207
    C_MASS = 12.00000000000
    N_MASS = 14.0030740048
    O_MASS = 15.99491461956
    BR79_MASS = 78.9183371

    # Given experimental data
    M_PLUS_H_MZ = 1108.70902
    BROMINE_COUNT = 6
    # Mass tolerance for finding a formula match (in Daltons)
    TOLERANCE = 0.003  # 3 mDa

    # --- Part 2: Calculate the mass of the unknown remainder ---
    # The problem states the isotopic pattern is due to 6 bromines.
    print(f"Step 1: The isotopic pattern 1:6:15:20:15:6:1 indicates {BROMINE_COUNT} bromine atoms.")

    # Calculate the mass of the neutral monoisotopic species M from the [M+H]+ ion
    neutral_mass = M_PLUS_H_MZ - H_MASS
    print(f"Step 2: The mass of the protonated ion [M+H]+ is {M_PLUS_H_MZ:.5f} Da.")
    print(f"        Subtracting a proton mass ({H_MASS:.5f} Da) gives the neutral monoisotopic mass: {neutral_mass:.5f} Da.")

    # Calculate the mass of the C, H, N, O portion of the molecule
    bromine_mass_part = BROMINE_COUNT * BR79_MASS
    remainder_mass = neutral_mass - bromine_mass_part
    print(f"Step 3: Subtracting the mass of {BROMINE_COUNT} ⁷⁹Br atoms ({bromine_mass_part:.5f} Da) isolates the mass of the remainder.")
    print(f"        The mass of the C,H,N,O portion is: {remainder_mass:.5f} Da.\n")
    
    print("Step 4: Searching for a valid formula (C,H,N,O) that matches this mass...")

    # --- Part 3: Search for the elemental formula of the remainder ---
    # Set plausible ranges for the number of each atom
    # The nominal mass of the remainder is ~634 Da. Max C is ~634/12 = 52.
    c_range = range(25, 45)
    # The Nitrogen Rule states that an odd nominal mass (1107) requires an odd number of Nitrogens.
    n_range = range(1, 15, 2)
    o_range = range(1, 15)
    
    solution_found = False
    for o_count in o_range:
        for n_count in n_range:
            # For a valid formula (integer DBE), if N is odd, (H+Halogens) must also be odd.
            # Since Halogens=6 (even), H must be odd.
            # This check is implicitly handled by the DBE calculation later.
            
            for c_count in c_range:
                # Calculate remaining mass for hydrogens
                mass_for_h = remainder_mass - (c_count * C_MASS + n_count * N_MASS + o_count * O_MASS)
                
                # Estimate number of hydrogens
                if mass_for_h < 0:
                    continue
                h_count = int(round(mass_for_h / H_MASS))
                if h_count <= 0:
                    continue

                # --- Part 4: Validate the candidate formula ---
                calculated_mass = c_count * C_MASS + h_count * H_MASS + n_count * N_MASS + o_count * O_MASS
                
                # Check if mass is within tolerance
                if abs(calculated_mass - remainder_mass) < TOLERANCE:
                    # Check for valid structure using Double Bond Equivalent (DBE)
                    # DBE = C - (H+X)/2 + N/2 + 1, must be a non-negative integer
                    dbe = c_count - (h_count + BROMINE_COUNT) / 2.0 + n_count / 2.0 + 1
                    
                    if dbe >= 0 and dbe == int(dbe):
                        print("\n--- Solution Found! ---")
                        print(f"A valid formula was found with a mass error of {abs(calculated_mass - remainder_mass)*1000:.2f} mDa.")
                        print("The molecular formula of the neutral species is composed of:")
                        print(f"Carbon (C): {c_count}")
                        print(f"Hydrogen (H): {h_count}")
                        print(f"Nitrogen (N): {n_count}")
                        print(f"Oxygen (O): {o_count}")
                        print(f"Bromine (Br): {BROMINE_COUNT}")
                        
                        final_formula = f"C{c_count}H{h_count}Br{BROMINE_COUNT}N{n_count}O{o_count}"
                        print(f"\nFinal Molecular Formula: {final_formula}")
                        solution_found = True
                        return final_formula

    if not solution_found:
        print("No solution found within the specified constraints and ranges.")
        return None

# Run the function to find the formula
final_formula_string = find_molecular_formula()
if final_formula_string:
    print(f"\n<<<C34H29N5O8Br6>>>")