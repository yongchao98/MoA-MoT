def solve():
    """
    This function determines for which k in {1, 2, ..., 9} the rational homotopy group
    pi_k(X) tensor Q vanishes, where X = S^4 wedge CP^2.

    The vanishing of the k-th rational homotopy group corresponds to the absence of
    a generator of degree k in the minimal Sullivan model of the space X.

    Based on the analysis:
    - The rational cohomology algebra is H*(X; Q) = Q[a,b]/<a^2, b^3, ab>.
    - Constructing the minimal model for this algebra reveals generators in the following degrees:
      - Degree 2: for the cohomology generator b.
      - Degree 4: for the cohomology generator a.
      - Degree 5: one to kill the cocycle corresponding to b^3=0, and another to kill the cocycle for ab=0.
      - Degree 7: to kill the cocycle corresponding to a^2=0.
    - Thus, the non-vanishing rational homotopy groups are for k = 2, 4, 5, 7.
    - We are looking for the k where the group *does* vanish.
    """

    all_k = set(range(1, 10))
    non_vanishing_k = {2, 4, 5, 7}
    vanishing_k = sorted(list(all_k - non_vanishing_k))
    
    # The final answer should be a comma-separated string of numbers.
    print(",".join(map(str, vanishing_k)))

solve()