def solve_electron_energy_problem():
    """
    This function calculates the minimum energy of electron 1 required for the
    transition of electron 2 from band II to band I, following the plan outlined above.
    """
    # We represent the bandgap energy E_g symbolically.
    Eg_str = "E_g"

    print("Step 1: State the conservation laws for the system.")
    print("Energy Conservation: E1_i + E2_i = E1_f + E2_f")
    print("Momentum Conservation: k1_i + k2_i = k1_f + k2_f\n")

    print("Step 2: Apply the condition for minimum required energy.")
    print("The transition is easiest when electron 2 starts at the top of band II, where k2_i = 0 and E2_i = 0.")
    print("The conservation laws become:")
    print("E1_i = E1_f + E2_f")
    print("k1_i = k1_f + k2_f\n")

    print("Step 3: Substitute the energy expressions and solve.")
    print("Using E_I = E_g + C*k^2, the energy equation is:")
    print(f"E_g + C*k1_i^2 = (E_g + C*k1_f^2) + (E_g + C*k2_f^2)")
    print(f"This simplifies to: C*k1_i^2 = {Eg_str} + C*(k1_f^2 + k2_f^2)\n")

    print("Step 4: Combine energy and momentum conservation.")
    print("Substitute k1_i = k1_f + k2_f into the equation from Step 3:")
    print("C*(k1_f + k2_f)^2 = E_g + C*(k1_f^2 + k2_f^2)")
    print("C*(k1_f^2 + k2_f^2 + 2*k1_f.k2_f) = E_g + C*(k1_f^2 + k2_f^2)")
    print("This gives the condition on the final states: 2*C*(k1_f.k2_f) = E_g\n")

    print("Step 5: Express the initial energy E1_i and minimize it.")
    print("From Step 3, we can write E1_i as:")
    print(f"E1_i = E_g + C*k1_i^2 = E_g + ({Eg_str} + C*(k1_f^2 + k2_f^2))")
    initial_energy_coeff = 2
    print(f"E1_i = {initial_energy_coeff} * {Eg_str} + C*(k1_f^2 + k2_f^2)\n")

    print("To minimize E1_i, we must minimize the term C*(k1_f^2 + k2_f^2) subject to the condition from Step 4.")
    print("This minimization occurs when the final electrons have equal momentum magnitudes (|k1_f| = |k2_f|) and move in the same direction.")
    print(f"In this case, the minimum value of C*(k1_f^2 + k2_f^2) is exactly E_g.")
    kinetic_part_coeff = 1
    print(f"So, (C*(k1_f^2 + k2_f^2))_min = {kinetic_part_coeff} * {Eg_str}\n")

    print("Step 6: Calculate the final result.")
    final_coeff = initial_energy_coeff + kinetic_part_coeff
    print("Substituting this minimum value back into the equation for E1_i:")
    print(f"E1_min = {initial_energy_coeff} * {Eg_str} + {kinetic_part_coeff} * {Eg_str}")
    print(f"E1_min = {final_coeff} * {Eg_str}")


solve_electron_energy_problem()