def find_upper_bound_relation():
    """
    This script derives and prints the upper bound for the maximum norm (k_{k,inf})
    in relation to the covolume (V) for a 2-dimensional lattice.
    """

    # The problem concerns lattices from real quadratic fields, so the dimension n is 2.
    n = 2
    V_symbol = "V"
    k_symbol = "k_{k,inf}"

    print("The upper bound is derived from Minkowski's theorem on convex bodies.")
    print(f"For a lattice of dimension 'n' and covolume '{V_symbol}', the maximum norm of the shortest non-zero vector (denoted '{k_symbol}') satisfies the inequality:")
    print(f"{k_symbol}^n <= {V_symbol}")
    print("")

    print(f"In the context of squarefree natural numbers (N), the associated lattices are 2-dimensional, so we set n = {n}.")
    print("Substituting n into the general inequality, we get:")
    print(f"{k_symbol}^{n} <= {V_symbol}")
    print("")

    print("To find the upper bound for the norm, we solve for it by taking the n-th root of both sides.")
    print(f"This results in: {k_symbol} <= {V_symbol}^(1/{n})")
    print("")

    print("So, the final equation expressing the upper bound is:")
    # Print each part of the final equation as requested
    print(k_symbol, end=' ')
    print("<=", end=' ')
    print(f"sqrt({V_symbol})")

find_upper_bound_relation()