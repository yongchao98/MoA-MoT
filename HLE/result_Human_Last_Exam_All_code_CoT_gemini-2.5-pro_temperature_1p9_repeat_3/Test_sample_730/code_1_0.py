def solve_parrot_problem():
    """
    This function formulates a calculation for a parrot to estimate the mass of a rock,
    adhering to the parrot's calculation constraints, and determines the final answer.
    """
    # The calculation is mass = density * (4/3) * pi * radius^3.
    # We will use fractional approximations for each term.

    # 1. Define the fractions based on the problem statement.
    # density = 0.9 = 9/10
    density_num = 9
    density_den = 10

    # volume constant = 4/3
    vol_const_num = 4
    vol_const_den = 3

    # radius = 0.5 = 1/2
    radius_num = 1
    radius_den = 2

    # 2. Approximate pi with a simple fraction.
    # As explained in the plan, pi ~ 3 is a valid approximation.
    # We represent it as 3/1.
    pi_approx_num = 3
    pi_approx_den = 1

    # 3. Assemble the list of all integers used in the calculation.
    all_integers = [
        density_num, density_den,
        vol_const_num, vol_const_den,
        pi_approx_num, pi_approx_den,
        radius_num, radius_den
    ]

    # Find the largest integer.
    largest_integer = max(all_integers)

    # 4. Print the calculation steps for the user as requested.
    # The prompt asks to output each number in the final equation.
    print("Yes, the parrot can estimate the mass with an error of at most 10%.")
    print("The calculation to perform is:")
    print(f"mass = ({density_num} / {density_den}) * ({vol_const_num} / {vol_const_den}) * {pi_approx_num} * ({radius_num} / {radius_den})^3")

    # For verification, calculate the result of this estimation.
    estimated_mass = (density_num / density_den) * (vol_const_num / vol_const_den) * pi_approx_num * (radius_num / radius_den)**3
    print(f"\nThe result of this calculation is {estimated_mass:.2f} kg.")
    
    # 5. Formulate and print the final answer.
    answer = f"Y{largest_integer}"
    print(f"The largest integer used in this calculation is {largest_integer}.")
    print(f"<<<{answer}>>>")

solve_parrot_problem()