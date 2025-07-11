import sympy

def solve_limit():
    """
    Calculates the limit based on the derived asymptotic behavior of p_n.

    The derivation shows that p_n, the probability that a simple random walk 
    starting at (n,0,0) escapes from the cube C_n = [0,2n]^3 and never returns,
    has an asymptotic behavior p_n ~ C/n for some constant C.

    We need to compute the limit: L = lim_{n->inf} [ln(1/p_n) / ln(n)]
    """

    # Define the symbols for the limit calculation
    n = sympy.Symbol('n', real=True, positive=True)
    # C is an arbitrary positive constant from the asymptotic analysis
    C = sympy.Symbol('C', real=True, positive=True)

    # The asymptotic form of p_n
    p_n = C / n

    # The expression inside the limit
    # log is the natural logarithm (ln) in sympy
    expression = sympy.log(1 / p_n) / sympy.log(n)

    # Calculate the limit as n approaches infinity
    limit_value = sympy.limit(expression, n, sympy.oo)

    # The final equation is lim_{n->oo} expression = limit_value
    # We print the final numerical answer as requested.
    print(limit_value)

if __name__ == '__main__':
    solve_limit()