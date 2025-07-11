import math

def solve_molecular_formula():
    """
    Determines the molecular formula of a neutral species from UPLC-MS data.
    """
    # --- Define constants for monoisotopic masses ---
    MASS_PROTONATED_ION = 1108.70902
    MASS_H = 1.007825
    MASS_C = 12.000000
    MASS_N = 14.003074
    MASS_O = 15.994915
    MASS_79Br = 78.918337
    
    # --- Step 1: Analyze the isotopic pattern ---
    print("Step 1: Determine the number of Bromine atoms from the isotopic pattern.")
    print("The isotopic distribution 1:6:15:20:15:6:1 corresponds to the binomial expansion of (a+b)^6.")
    print("Given that Bromine's isotopes (79Br and 81Br) have a ~1:1 abundance and differ by 2 amu, this pattern indicates the presence of 6 Bromine atoms.")
    num_br = 6
    print(f"Conclusion: The molecule contains {num_br} Bromine atoms.\n")

    # --- Step 2: Calculate the mass of the neutral molecule [M] ---
    print("Step 2: Calculate the mass of the neutral molecule [M].")
    mass_neutral = MASS_PROTONATED_ION - MASS_H
    print(f"The observed m/z for the protonated ion [M+H]+ is {MASS_PROTONATED_ION}.")
    print(f"Subtracting the mass of a proton (H = {MASS_H} Da), we get the mass of the neutral molecule [M].")
    print(f"Mass [M] = {MASS_PROTONATED_ION} - {MASS_H} = {mass_neutral:.6f} Da\n")

    # --- Step 3: Calculate the mass of the non-bromine remainder ---
    print("Step 3: Calculate the mass of the remainder (non-bromine part).")
    mass_br_total = num_br * MASS_79Br
    mass_remainder = mass_neutral - mass_br_total
    print(f"The lowest m/z peak corresponds to the monoisotopic mass, containing only the lightest isotopes (e.g., 79Br).")
    print(f"Mass of six 79Br atoms = {num_br} * {MASS_79Br} = {mass_br_total:.6f} Da.")
    print(f"Mass of Remainder = Mass [M] - Mass [6 x Br] = {mass_neutral:.6f} - {mass_br_total:.6f} = {mass_remainder:.6f} Da.\n")

    # --- Step 4: Determine the formula of the remainder (C, H, N, O part) ---
    print("Step 4: Determine the molecular formula for the remainder.")
    print(f"The nominal mass of the neutral molecule is {math.floor(mass_neutral)}, which is odd.")
    print("According to the Nitrogen Rule, an odd nominal mass implies an odd number of nitrogen atoms.")
    print("Searching for a combination of C, H, N, O that matches the remainder mass of {:.6f} Da...".format(mass_remainder))
    
    # Search for the formula of the remainder within plausible ranges
    target_mass = mass_remainder
    tolerance = 0.003  # 3 mDa tolerance
    solution_found = False
    
    for c in range(30, 40):
        if solution_found: break
        for o in range(5, 12):
            if solution_found: break
            # Iterate through odd numbers of N according to the Nitrogen Rule
            for n in range(1, 11, 2):
                mass_cno = c * MASS_C + o * MASS_O + n * MASS_N
                # Calculate the number of hydrogens needed
                h_mass_needed = target_mass - mass_cno
                num_h = round(h_mass_needed / MASS_H)

                if num_h <= 0:
                    continue

                calculated_mass = mass_cno + num_h * MASS_H
                
                # Check if the calculated mass is within tolerance
                if abs(calculated_mass - target_mass) < tolerance:
                    # Check for valid degree of unsaturation (must be integer or half-integer)
                    if (num_h - n) % 2 == 0:
                        dbe = c - num_h / 2 + n / 2 + 1
                        if dbe >= 0:
                            print(f"\nFound a plausible match for the remainder: C{c}H{num_h}N{n}O{o}")
                            print(f"  - Calculated Mass: {calculated_mass:.6f} Da")
                            print(f"  - Mass Error: {abs(calculated_mass - target_mass)*1e6/target_mass:.2f} ppm")
                            
                            # --- Step 5: Final Answer ---
                            print("\nStep 5: Assemble the full molecular formula.")
                            print("Combining the bromine atoms with the remainder gives the final formula.")
                            print("\nThe molecular formula of the neutral species is:")
                            print(f"C{c}H{num_h}Br{num_br}N{n}O{o}")
                            final_formula = f"C{c}H{num_h}Br{num_br}N{n}O{o}"
                            
                            # Final verification
                            total_mass_check = calculated_mass + mass_br_total
                            total_protonated_check = total_mass_check + MASS_H
                            print(f"\nVerification:")
                            print(f"Calculated [M+H]+ Mass: {total_protonated_check:.6f} Da")
                            print(f"Observed [M+H]+ Mass:   {MASS_PROTONATED_ION:.6f} Da")
                            print(f"Error: {abs(total_protonated_check - MASS_PROTONATED_ION)*1e6/MASS_PROTONATED_ION:.2f} ppm")

                            print(f"\n<<<C{c}H{num_h}Br{num_br}N{n}O{o}>>>")
                            solution_found = True
                            break
    if not solution_found:
        print("No plausible formula was found within the specified search range.")

# Run the solver
solve_molecular_formula()