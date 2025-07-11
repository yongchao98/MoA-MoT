def solve_markov_chain_problem():
    """
    This function calculates and prints the symbolic expression for the supremum of alpha.
    The derivation is based on the diffusion approximation of the Markov chain.

    The final result depends on the parameter 'c' from the problem description.
    """

    # These are the numerical constants from the derived formula for the supremum of alpha.
    numerator_constant = 1
    numerator_c_coefficient = 4
    denominator_val = 2

    # The problem asks for sup{alpha: E[tau^alpha] < infinity}.
    # Our derived formula is (1 + 4*c) / 2.
    # The code below prints this formula in a structured way.

    print("The final equation for the supremum of alpha is derived from the diffusion approximation:")
    print(f"sup(alpha) = ({numerator_constant} + {numerator_c_coefficient}*c) / {denominator_val}")


solve_markov_chain_problem()
