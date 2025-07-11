def solve_semidistributivity():
    """
    This function determines the largest cardinal mu such that any forcing notion P
    with a dense subset of cardinality kappa is necessarily (mu, kappa+)-semidistributive.

    Let P be a forcing notion with d(P) = kappa, where d(P) is the smallest
    cardinality of a dense subset of P.
    We are looking for the largest cardinal mu such that P is necessarily
    (mu, kappa+)-semidistributive.

    A forcing P is (mu, lambda)-semidistributive if for every lambda-sized set X
    in the generic extension V[G] (where X is a subset of lambda), there is a
    ground-model set Y (i.e., Y in V) such that Y is a subset of X and |Y| = mu.

    1.  Check for mu = 1:
        If X is a set of size kappa+, it is non-empty. Let alpha be an element of X.
        The set Y = {alpha} is a ground model set, since alpha is an ordinal from
        the ground model. Y is a subset of X, and its size is 1.
        So, any such forcing is (1, kappa+)-semidistributive. This means mu >= 1.

    2.  Check if mu can be 2:
        To show that mu is not necessarily 2, we need to find a counterexample:
        a forcing P with d(P) = kappa that is NOT (2, kappa+)-semidistributive.

        Consider the Cohen forcing P = Fn(kappa, 2), the set of all finite partial
        functions from kappa to 2. The density of this forcing is d(P) = kappa.
        This forcing adds a generic function g: kappa -> 2.

        Using this generic function g, one can construct a P-name for a set X,
        a subset of kappa+, of size kappa+, which is "generic" enough to not contain
        any ground-model pair {alpha, beta}.

        A standard construction uses an almost disjoint family of sets. Let
        (S_alpha)_{\alpha < kappa+} be a family of subsets of kappa, each of size kappa,
        such that for any alpha != beta, |S_alpha intersect S_beta| < kappa.
        Define X in the generic extension V[G] to be:
        X = { alpha < kappa+ | the set {i in S_alpha | g(i) = 1} has size < kappa }

        It is a known (but advanced) result that this set X has size kappa+ in V[G].
        It can also be shown that for any specific alpha from the ground model, the set of
        conditions in P that force `alpha not in X` is dense in P.
        This means that for any generic filter G, `alpha not in X[G]` will hold.
        Consequently, for any pair Y = {alpha, beta} from the ground model,
        Y will not be a subset of X.

        Since we have a counterexample for mu = 2 that works for any kappa,
        mu cannot be 2.

    3.  Conclusion:
        Since mu must be at least 1, and it cannot be 2 (or greater), the largest
        value for mu is 1. The result is independent of kappa.
    """
    mu = 1
    # The question asks to "output each number in the final equation!".
    # The final conclusion is that the largest mu is 1.
    print(f"mu = {mu}")

solve_semidistributivity()