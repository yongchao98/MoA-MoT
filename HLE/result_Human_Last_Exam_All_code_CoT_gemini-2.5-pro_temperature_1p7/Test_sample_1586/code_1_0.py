def solve_markov_problem():
    """
    This function prints the derived formula for the supremum of alpha.
    The problem is solved analytically, and this code serves to present the result.
    The final answer is a formula dependent on the constant 'c'.
    """

    # The numbers present in the final derived formula.
    number_one = 1
    number_four = 4
    number_two = 2

    # The final equation for the supremum of alpha is sup(alpha) = (1 + 4*c) / 2.
    # We construct and print this equation string by string to highlight each component.
    
    print("The analytical solution for the supremum of alpha is given by the following equation:")
    # We use an f-string to embed the numbers into the output string.
    # The variable 'c' is treated as a symbolic placeholder.
    print(f"sup(alpha) = ({number_one} + {number_four} * c) / {number_two}")

solve_markov_problem()