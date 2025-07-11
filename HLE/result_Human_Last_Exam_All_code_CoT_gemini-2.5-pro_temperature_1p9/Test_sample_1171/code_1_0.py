def solve_homotopy_groups():
    """
    This function determines for which k in {1, ..., 9} the rational homotopy group
    pi_k(S^4 v CP^2) vanishes.
    """
    
    # K is the set of integers we are interested in.
    K = set(range(1, 10))

    # Degrees of generators for the minimal Sullivan model of CP^2
    gen_CP2 = {2, 5}
    
    # Degrees of generators for the minimal Sullivan model of S^4
    gen_S4 = {4, 7}

    # Degrees of linking generators
    # H_tilde(CP^2) has generators a (deg 2) and a^2 (deg 4)
    # H_tilde(S^4) has generator b (deg 4)
    # Product a*b has degree 2+4=6, needs a linking generator of degree 6-1=5.
    # Product a^2*b has degree 4+4=8, needs a linking generator of degree 8-1=7.
    gen_linking = {5, 7}

    # The set of k for which the rational homotopy group is non-zero
    # is the union of the degrees of all generators.
    non_vanishing_k = gen_CP2.union(gen_S4).union(gen_linking)

    # The set of k for which the group vanishes is the complement in K.
    vanishing_k = sorted(list(K - non_vanishing_k))

    # Print the result in the specified format.
    print(','.join(map(str, vanishing_k)))

solve_homotopy_groups()