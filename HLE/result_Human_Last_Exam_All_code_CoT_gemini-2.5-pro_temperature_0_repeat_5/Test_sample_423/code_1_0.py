def solve_structure():
    """
    Performs calculations to help determine the structure of substance X.
    1. Calculates the C:H ratio from combustion data.
    2. Identifies possible molecular formulas based on the molar mass range.
    3. Prints the derived structures of intermediates A and B.
    4. Proposes a final structure for X based on the available clues.
    """
    # Molar masses
    M_C = 12.011
    M_H = 1.008
    M_O = 15.999
    M_CO2 = 44.01
    M_H2O = 18.015

    # Combustion data
    mass_co2 = 0.7472
    mass_h2o = 0.1834

    # Molar mass estimation
    molar_mass_x_approx = 150
    error_margin = 0.10
    molar_mass_min = molar_mass_x_approx * (1 - error_margin)
    molar_mass_max = molar_mass_x_approx * (1 + error_margin)

    # Step 1: Calculate moles of C and H
    moles_c = mass_co2 / M_CO2
    moles_h = (mass_h2o / M_H2O) * 2

    # Step 2: Find the C:H ratio
    # The ratio is moles_c : moles_h
    # To find the simplest integer ratio, divide by the smaller number
    # and then multiply to get integers.
    # H/C ratio = moles_h / moles_c = 2.036e-2 / 1.698e-2 = 1.199 ~= 1.2 = 6/5
    c_ratio = 5
    h_ratio = 6

    print("Step 1: Analysis of X")
    print(f"Moles of C = {moles_c:.4f} mol")
    print(f"Moles of H = {moles_h:.4f} mol")
    print(f"The simplest integer ratio C:H is {c_ratio}:{h_ratio}\n")

    # Step 3: Find possible molecular formulas
    print("Possible molecular formulas for X (C, H, O) based on Molar Mass ~150 (135-165 g/mol):")
    possible_formulas = []
    # Check for multiples of the empirical C:H unit (C5H6)
    for n in range(1, 3): # Multiplier for C5H6
        for o_atoms in range(1, 10): # Number of oxygen atoms
            c_atoms = c_ratio * n
            h_atoms = h_ratio * n
            formula = f"C{c_atoms}H{h_atoms}O{o_atoms}"
            mass = c_atoms * M_C + h_atoms * M_H + o_atoms * M_O
            if molar_mass_min <= mass <= molar_mass_max:
                possible_formulas.append((formula, mass))
                print(f"- Formula: {formula}, Molar Mass: {mass:.2f} g/mol")
    
    print("\nStep 2: Structure of B")
    # Based on elemental analysis (C:0.5, H:0.1, O:0.4 -> C5H12O3) and reaction with HBr
    structure_b = "Pentane-1,3,5-triol"
    print(f"The structure of B is {structure_b}.\n")

    print("Step 3: Structure of A")
    # Based on reduction to B
    structure_a = "Pentane-1,5-dial-3-one"
    print(f"The structure of A is {structure_a}.\n")

    print("Step 4: Proposed Structure of X")
    print("The provided data contains contradictions, particularly between the combustion analysis/molar mass and the ozonolysis reaction.")
    print("The ozonolysis of X to a single C5 product (A) from a C10 starting material (X) implies X is a symmetrical dimer.")
    print("This path leads to a molecular formula of C10H12O4 (M=196), which contradicts the estimated molar mass.")
    print("Another path suggests X is a resorcinol (1,3-dihydroxybenzene) derivative, as its ozonolysis can yield A.")
    print("However, no resorcinol derivative fits all the data (formula C10H12O2, Tollen's positive).")
    print("Given the strong evidence for a phenolic aldehyde structure that is also a resorcinol derivative, a plausible, albeit imperfect, answer is a known compound that fits these chemical properties, even if the formula doesn't perfectly match the flawed quantitative data.")
    
    final_structure_x = "2,4-dihydroxy-6-methylbenzaldehyde (Orsellinic aldehyde)"
    print(f"\nA plausible structure for X, assuming errors in the quantitative data, is {final_structure_x}.")
    print("This structure is a resorcinol derivative and an aldehyde, fitting the key chemical tests.")

solve_structure()
<<<2,4-dihydroxy-6-methylbenzaldehyde>>>