def solve_markov_chain_moment():
    """
    This function prints the derived supremum for alpha.
    The problem asks to find sup{alpha: E[tau^alpha] < infinity} for a given Markov chain.
    
    The final equation derived from the analysis is:
    sup_alpha = 4 * c
    
    where c is the positive constant from the problem description.
    """
    # The equation is sup_alpha = 4 * c.
    # The number in the equation is 4.
    number_in_equation = 4
    
    # We print the final result as an expression in terms of c.
    # The prompt requests to "output each number in the final equation",
    # so we construct the output string to show the number clearly.
    print(f"The supremum of alpha is: {number_in_equation}*c")

solve_markov_chain_moment()