def solve():
    """
    Finds the values of k in {1, 2, ..., 9} for which the rational homotopy group
    pi_k(S^4 v CP^2) vanishes.
    Based on the analysis of the minimal Sullivan model for the cohomology algebra of X.
    d_k = dim(pi_k(X) tensor Q).
    """
    # Dimensions d_k = dim(pi_k(X) tensor Q) for k=1 to 9
    # These are derived from the construction of the minimal model.
    # d_k is the number of generators in degree k.
    d = {
        1: 0,  # X is path-connected and simply connected
        2: 1,  # Corresponds to H^2
        3: 0,  # No generator needed
        4: 1,  # Corresponds to the second generator of H^4
        5: 2,  # To kill cocycles in degree 6 (v_2^3 and v_2*v_4)
        6: 0,  # No new generator needed
        7: 1,  # To kill cocycle in degree 8 (v_4^2)
        8: 0,  # No new generator needed
        9: 1,  # To kill cocycle in degree 10 (v_2^5)
    }

    vanishing_k = []
    for k in range(1, 10):
        if d.get(k) == 0:
            vanishing_k.append(str(k))

    print(",".join(vanishing_k))

solve()