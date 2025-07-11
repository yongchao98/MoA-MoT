import math

def solve_mass_ratio():
    """
    Calculates and explains the leading-order asymptotic mass ratio
    in the CP(N-1) model.
    """
    print("Step 1: State the mass spectrum formula for k-particle bound states in the CP(N-1) model.")
    print("The formula is: M_k = m * sin(pi * k / N) / sin(pi / N)\n")

    print("Step 2: Identify the masses of the two excitations in question.")
    print("The lightest excitation (k=1) has mass: M_1 = m * sin(pi / N) / sin(pi / N) = m")
    print("The subsequent higher excitation (k=2) has mass: M_2 = m * sin(2 * pi / N) / sin(pi / N)\n")

    print("Step 3: Formulate the mass ratio R = M_2 / M_1.")
    print("R = [m * sin(2 * pi / N) / sin(pi / N)] / m")
    print("R = sin(2 * pi / N) / sin(pi / N)\n")

    print("Step 4: Simplify the ratio using the identity sin(2x) = 2*sin(x)*cos(x).")
    print("Let x = pi / N. The ratio becomes:")
    print("R = [2 * sin(pi / N) * cos(pi / N)] / sin(pi / N)")
    print("R = 2 * cos(pi / N)\n")

    print("Step 5: Determine the asymptotic limit as N approaches infinity.")
    print("As N -> infinity, the term (pi / N) -> 0.")
    print("Therefore, cos(pi / N) -> cos(0), which is 1.\n")
    
    print("Step 6: State the final result for the asymptotic mass ratio.")
    # The final equation with each number
    part_1 = 2
    part_2 = int(math.cos(0))
    final_result = part_1 * part_2
    print(f"The final ratio is: {part_1} * {part_2} = {final_result}")

solve_mass_ratio()
<<<2>>>