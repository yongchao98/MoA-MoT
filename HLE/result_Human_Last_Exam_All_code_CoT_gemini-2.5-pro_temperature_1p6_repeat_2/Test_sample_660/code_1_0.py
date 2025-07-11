def solve():
    """
    This function derives and prints the condition under which companies with more
    truthful accounts choose more lenient auditors.

    The cost for a firm of type theta choosing an auditor of type x is:
    C(theta, 0) = t(theta)
    C(theta, 1) = theta*t(theta) + (1-theta)*p(1)*(t(theta)+F(theta))

    The condition "more truthful firms choose more lenient auditors" implies a
    sorting where x*(1) <= x*(0). This holds if the incremental cost of
    hiring a strict auditor, Delta_C(theta) = C(theta, 1) - C(theta, 0),
    is non-decreasing in theta. For theta in {0, 1}, this means
    Delta_C(1) >= Delta_C(0).

    Delta_C(theta) = (1-theta) * [p(1)*(t(theta) + F(theta)) - t(theta)]

    Delta_C(1) = 0
    Delta_C(0) = p(1)*(t(0) + F(0)) - t(0)

    The condition becomes:
    0 >= p(1)*(t(0) + F(0)) - t(0)
    which simplifies to:
    t(0) >= p(1)*(t(0) + F(0))
    """
    
    # We will print the final condition, term by term.
    # The terms are symbolic representations of the functions at specific points.
    t_0 = "t(0)"
    p_1 = "p(1)"
    F_0 = "F(0)"
    
    # Printing each component of the final equation: t(0) >= p(1) * (t(0) + F(0))
    print(t_0, ">=", p_1, "*", "(", t_0, "+", F_0, ")")

solve()