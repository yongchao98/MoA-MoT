def solve_homotopy_vanishing():
    """
    Calculates for which k in {1, ..., 9} the rational homotopy group
    pi_k(S^4 v CP^2) x Q vanishes.
    """
    
    # The set of k values we are interested in.
    k_values = set(range(1, 10))
    
    # For S^4, the non-zero rational homotopy groups for k<=9 are for k=4 and k=7.
    # pi_k(S^4) x Q is non-zero for k in {4, 7}.
    non_vanishing_s4 = {4, 7}
    
    # For CP^2, the non-zero rational homotopy groups are pi_2 and pi_5.
    # pi_k(CP^2) x Q is non-zero for k in {2, 5}.
    non_vanishing_cp2 = {2, 5}
    
    # For the wedge sum X = S^4 v CP^2, pi_k(X) x Q is the direct sum of the
    # rational homotopy groups of the components. It is non-zero if at least
    # one of the components has a non-zero group.
    non_vanishing_x = non_vanishing_s4.union(non_vanishing_cp2)
    
    # The group vanishes for k values that are not in the non_vanishing_x set.
    vanishing_k = k_values.difference(non_vanishing_x)
    
    # Sort the final list for a clean output.
    result = sorted(list(vanishing_k))
    
    # Print the explanation of the "equation"
    print(f"The set of k in {{1, ..., 9}} is {sorted(list(k_values))}.")
    print(f"pi_k(S^4) x Q is non-zero for k in {sorted(list(non_vanishing_s4))}.")
    print(f"pi_k(CP^2) x Q is non-zero for k in {sorted(list(non_vanishing_cp2))}.")
    print(f"Thus, pi_k(S^4 v CP^2) x Q is non-zero for k in the union, {sorted(list(non_vanishing_x))}.")
    print(f"The group vanishes for k in the complement of this set.")
    
    # Print the final answer in the required format.
    print(",".join(map(str, result)))

solve_homotopy_vanishing()