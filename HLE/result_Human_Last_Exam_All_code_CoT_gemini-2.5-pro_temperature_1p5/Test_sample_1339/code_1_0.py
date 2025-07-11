def solve_group_theory_problem():
    """
    This function formulates and prints the solution to the given problem.
    """
    # (a) Does there exist a unique minimal group G_hat?
    # Based on the theory of p-completion of groups, such a group exists, is unique
    # up to isomorphism over G, and is minimal.
    answer_a = "Yes"

    # (b) What is the maximum possible derived length of G_hat?
    # As determined by the step-by-step reasoning, the maximum possible derived length
    # of G_hat is n, where n is the number of abelian factors in the subnormal series.
    # The final equation for the maximum derived length is `d(G_hat) = n`.
    answer_b = "n"

    print(f"(a) {answer_a}; (b) {answer_b}")

solve_group_theory_problem()