def solve_mass_problem():
    """
    This function calculates the total mass based on the derived formula.
    The problem asks for the total mass of (q_v * (q - 1)) / (q_v - 1) * mu.
    Given q_v = q, this simplifies to q * mu.

    Our derivation shows that in the stable case (n -> infinity), the mass mu is 1/q.
    Therefore, the total mass is q * (1/q) = 1.

    To demonstrate this and output the numbers in the final equation as requested,
    we use an example value for q. The result is independent of this choice.
    """
    # Let's choose an example value for q, for instance 2.
    # The characteristic q must be a prime power > 1.
    q = 2

    # In the stable limit n -> infinity, the mass mu(X) is 1/q.
    mu_X = 1 / q

    # The quantity to compute is q * mu(X).
    total_mass = q * mu_X

    # Output the numbers in the final equation.
    print(f"The calculation for an example value q={q} is:")
    print(f"{q} * {mu_X} = {total_mass}")
    print("\nThe final result is independent of q and is equal to 1.")

solve_mass_problem()