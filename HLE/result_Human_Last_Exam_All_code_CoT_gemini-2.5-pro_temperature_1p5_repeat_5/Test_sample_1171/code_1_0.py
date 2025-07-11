def solve_homotopy_vanishing():
    """
    Solves for which k in {1, ..., 9} the rational homotopy group pi_k(X) tensor Q vanishes,
    where X = S^4 v CP^2.

    The rational homotopy groups of a space are non-zero only at degrees k corresponding
    to the degrees of the generators of its minimal Sullivan model.

    1. The rational cohomology of X is H*(X; Q) = Q[a,b]/(a^3, b^2, ab),
       with deg(a)=2, deg(b)=4.

    2. From this, we construct the minimal Sullivan model for X. The generators are found to have degrees:
       - One generator of degree 2 (for H^2)
       - One generator of degree 4 (for H^4)
       - Two generators of degree 5 (to kill the cycles a^3 and ab in cohomology)
       - One generator of degree 7 (to kill the cycle b^2 in cohomology)

    3. Thus, the non-vanishing rational homotopy groups pi_k(X) tensor Q are for k in {2, 4, 5, 7}.

    4. We are asked for which k in {1, ..., 9} the group *vanishes*. We take the complement
       of {2, 4, 5, 7} in the set {1, 2, 3, 4, 5, 6, 7, 8, 9}.
    """
    
    all_k = set(range(1, 10))
    non_vanishing_k = {2, 4, 5, 7}
    vanishing_k = sorted(list(all_k - non_vanishing_k))
    
    # Print the result in the required format
    print(','.join(map(str, vanishing_k)))

solve_homotopy_vanishing()