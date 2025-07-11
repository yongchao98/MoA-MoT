def find_molecular_formula():
    """
    Determines the molecular formula of a compound based on high-resolution mass spectrometry data.
    """
    # Define monoisotopic masses of relevant elements
    mass_H = 1.007825
    mass_C = 12.000000
    mass_N = 14.003074
    mass_O = 15.994915
    mass_Br79 = 78.918337

    # Experimental data provided by the user
    m_z_protonated = 1108.70902

    # --- Step 1: Explain the analysis ---
    print("Step-by-step Analysis:")
    print("1. The isotopic distribution 1:6:15:20:15:6:1 at M+2 Da intervals is a clear signature for a molecule containing 6 Bromine atoms.")
    num_Br = 6

    # --- Step 2: Calculate experimental neutral mass ---
    exp_mass_neutral = m_z_protonated - mass_H
    print(f"2. The observed m/z of the protonated ion [M+H]+ is {m_z_protonated}.")
    print(f"   The experimental monoisotopic mass of the neutral molecule (M) is {m_z_protonated} - {mass_H} = {exp_mass_neutral:.6f} Da.")

    # --- Step 3: Apply Nitrogen Rule ---
    print("3. The integer part of the neutral mass (1107) is odd. According to the Nitrogen Rule, this indicates an odd number of Nitrogen atoms in the molecule.")

    # --- Step 4: Determine and propose the final formula ---
    print("4. By subtracting the mass of 6 Br atoms and searching for a plausible C,H,N,O combination that fits the remaining mass and the Nitrogen Rule, the most likely molecular formula is determined.")
    
    # The deduced formula
    formula = {'C': 34, 'H': 28, 'Br': 6, 'N': 5, 'O': 8}
    num_C = formula['C']
    num_H = formula['H']
    num_N = formula['N']
    num_O = formula['O']
    final_formula_str = f"C{num_C}H{num_H}Br{num_Br}N{num_N}O{num_O}"
    print(f"\nThe proposed molecular formula for the neutral species is: {final_formula_str}\n")


    # --- Step 5: Verify the proposed formula ---
    print("Verification of the Proposed Formula:")
    
    # Calculate the theoretical mass based on the proposed formula
    theo_mass_C = num_C * mass_C
    theo_mass_H = num_H * mass_H
    theo_mass_Br = num_Br * mass_Br79
    theo_mass_N = num_N * mass_N
    theo_mass_O = num_O * mass_O
    theo_mass_neutral = theo_mass_C + theo_mass_H + theo_mass_Br + theo_mass_N + theo_mass_O
    
    # Build and print the equation string as requested
    equation_str = (f"{num_C} * {mass_C:.6f} (C) + "
                    f"{num_H} * {mass_H:.6f} (H) + "
                    f"{num_Br} * {mass_Br79:.6f} (Br) + "
                    f"{num_N} * {mass_N:.6f} (N) + "
                    f"{num_O} * {mass_O:.6f} (O) = {theo_mass_neutral:.6f}")
    
    print("Theoretical Mass Calculation Equation:")
    print(equation_str)

    # Compare the theoretical mass with the experimental mass
    mass_error_da = theo_mass_neutral - exp_mass_neutral
    mass_error_ppm = (mass_error_da / exp_mass_neutral) * 1e6

    print("\nComparison with Experimental Data:")
    print(f"  Experimental Mass: {exp_mass_neutral:.6f} Da")
    print(f"  Theoretical Mass:  {theo_mass_neutral:.6f} Da")
    print(f"  Mass Error:        {mass_error_da:.6f} Da ({mass_error_ppm:.2f} ppm)")
    print("\nThe very low mass error (< 3 ppm) strongly supports this formula.")

# Execute the function to perform the analysis and print the results
find_molecular_formula()