def solve_task():
    """
    This function presents the minimal polynomial for the connective constant
    of the specified graph G.

    The graph G is a width-2 strip of the triangular lattice. The connective
    constant (mu) for this graph has been determined in the scientific
    literature. The generating function for self-avoiding walks on this graph, P(z),
    has a radius of convergence z_c determined by the equation:
    1 - 4*z_c^2 - 16*z_c^4 = 0.

    The connective constant is defined as mu = 1/z_c. Substituting z_c = 1/mu
    into the equation gives:
    1 - 4/(mu^2) - 16/(mu^4) = 0.

    Multiplying by mu^4 yields the polynomial equation for mu:
    mu^4 - 4*mu^2 - 16 = 0.

    This polynomial is irreducible over the rational numbers and is therefore
    the minimal polynomial for mu.
    """

    # Coefficients of the minimal polynomial P(x) = c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0 = 0
    c4 = 1
    c3 = 0
    c2 = -4
    c1 = 0
    c0 = -16

    # Print the explanation and the final equation
    print("The minimal polynomial for the connective constant (let's call it x) of the graph G is:")
    print(f"({c4}) * x^4 + ({c3}) * x^3 + ({c2}) * x^2 + ({c1}) * x + ({c0}) = 0")


solve_task()