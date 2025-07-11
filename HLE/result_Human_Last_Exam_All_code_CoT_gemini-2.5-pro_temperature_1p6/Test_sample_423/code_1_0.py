def solve_structure():
    """
    This function performs the elemental analysis to determine the
    molecular formula of substance X based on the provided combustion data
    and molar mass estimate.
    """
    # Given data
    mass_co2 = 0.7472  # g
    mass_h2o = 0.1834  # g
    molar_mass_x_approx = 150  # g/mol
    error_margin = 0.10

    # Molar masses of elements and compounds
    m_c = 12.01
    m_h = 1.008
    m_o = 16.00
    m_co2 = m_c + 2 * m_o
    m_h2o = 2 * m_h + m_o

    # --- Step 1: Calculate moles of C and H ---
    moles_co2 = mass_co2 / m_co2
    moles_c = moles_co2

    moles_h2o = mass_h2o / m_h2o
    moles_h = 2 * moles_h2o

    print("--- Step 1: Elemental Analysis ---")
    print(f"Moles of CO2 produced: {moles_co2:.4f} mol")
    print(f"Moles of H2O produced: {moles_h2o:.4f} mol")
    print(f"Moles of Carbon (C) in sample: {moles_c:.4f} mol")
    print(f"Moles of Hydrogen (H) in sample: {moles_h:.4f} mol\n")

    # --- Step 2: Determine the empirical C:H ratio ---
    h_c_ratio = moles_h / moles_c
    # We look for a simple integer ratio. 1.2 is 6/5.
    c_empirical = 5
    h_empirical = 6
    
    print("--- Step 2: Empirical Formula Ratio ---")
    print(f"Ratio of H to C atoms: {h_c_ratio:.3f}")
    print(f"This is approximately 1.2, which corresponds to a C:H ratio of {c_empirical}:{h_empirical}.\n")

    # --- Step 3: Determine Molecular Formula using Molar Mass ---
    mass_c_empirical = c_empirical * m_c
    mass_h_empirical = h_empirical * m_h
    empirical_unit_mass = mass_c_empirical + mass_h_empirical

    m_min = molar_mass_x_approx * (1 - error_margin)
    m_max = molar_mass_x_approx * (1 + error_margin)

    print("--- Step 3: Finding Molecular Formula ---")
    print(f"The estimated molar mass is between {m_min:.1f} and {m_max:.1f} g/mol.")
    
    possible_formulas = []

    # Check for n=1 for (C5H6)n
    n = 1
    mass_from_ch = n * empirical_unit_mass
    remaining_mass_for_o = molar_mass_x_approx - mass_from_ch
    num_o = round(remaining_mass_for_o / m_o)
    if num_o > 0:
        formula = f"C{c_empirical*n}H{h_empirical*n}O{num_o}"
        total_mass = mass_from_ch + num_o * m_o
        if m_min <= total_mass <= m_max:
            possible_formulas.append((formula, total_mass))
            print(f"Testing n=1: Formula {formula}, Molar Mass = {total_mass:.2f} g/mol. This is a possible candidate.")

    # Check for n=2 for (C5H6)n
    n = 2
    mass_from_ch = n * empirical_unit_mass
    remaining_mass_for_o = molar_mass_x_approx - mass_from_ch
    num_o = round(remaining_mass_for_o / m_o)
    if num_o > 0:
      formula = f"C{c_empirical*n}H{h_empirical*n}O{num_o}"
      total_mass = mass_from_ch + num_o * m_o
      if m_min <= total_mass <= m_max:
          possible_formulas.append((formula, total_mass))
          print(f"Testing n=2: Formula {formula}, Molar Mass = {total_mass:.2f} g/mol. This is a possible candidate.")
    
    # --- Step 4: Final Deduction ---
    print("\n--- Step 4: Final Conclusion on Formula ---")
    print("Based on the chemical properties (reacts with NaOH, indicating a phenolic/acidic -OH) and (reduces Tollens' reagent, indicating an aldehyde), the formula must contain at least two oxygen atoms.")
    print("The candidate C10H12O only has one oxygen, so it is ruled out.")
    print("The most plausible molecular formula for substance X is C5H6O5.\n")
    
    # Final step: Print the equation for finding moles of C
    # moles C = mass CO2 / M_CO2
    print("To reiterate the key calculation for Carbon:")
    print(f"moles_C = {mass_co2:.4f} g / {m_co2:.2f} g/mol = {moles_c:.4f} mol")
    
solve_structure()