import sympy

def solve_hamiltonicity_threshold():
    """
    Calculates and prints the d-threshold for Hamiltonicity based on boosting
    the minimum degree to satisfy Dirac's theorem.

    The plan is as follows:
    1.  Start with a graph H_n with minimum degree δ(H_n) >= n/2 - η.
    2.  Use Dirac's theorem, which states that a graph G with n vertices is
        Hamiltonian if its minimum degree δ(G) is at least n/2.
    3.  Find the probability p such that the union graph H_n U G(n, p)
        satisfies Dirac's condition a.a.s.
    4.  The minimum degree of the union graph must be at least n/2.
        For a vertex v with the lowest possible degree in H_n (which is n/2 - η),
        its degree in the random graph G(n,p) must be at least η.
        So, we require δ(G(n, p)) >= η.
    5.  The known threshold for the minimum degree of a random graph G(n, p) to
        be at least k is p = (log(n) + (k-1)*log(log(n))) / n.
    6.  Substituting k = η gives the required probability threshold p.
    """
    
    # Define n and eta as symbolic variables for the formula.
    n = sympy.Symbol('n')
    eta = sympy.Symbol('eta')

    # The formula for the threshold is derived from the threshold for the
    # minimum degree of a random graph G(n, p) being at least eta.
    # The term '-1' is part of this established formula.
    p_threshold = (sympy.log(n) + (eta - 1) * sympy.log(sympy.log(n))) / n
    
    print("The d-threshold for Hamiltonicity is given by the formula for p:")
    
    # Create and pretty-print the equation p = ...
    # This fulfills the requirement to output each number in the final equation (here, the number 1).
    equation = sympy.Eq(sympy.Symbol('p'), p_threshold)
    sympy.pprint(equation, use_unicode=True)

solve_hamiltonicity_threshold()