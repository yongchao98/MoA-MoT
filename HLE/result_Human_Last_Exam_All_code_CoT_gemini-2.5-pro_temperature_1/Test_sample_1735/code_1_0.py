def solve():
    """
    This function prints the sorted letter choices of the true statements.
    Based on the analysis:
    A) False. The problem is NP-hard.
    B) False. Optimal solutions often interleave positive and negative numbers.
    C) True. The objective function |S_k| is symmetric with respect to negation.
    D) False. Verification is as hard as solving an NP-hard problem.
    E) False. Non-uniqueness can occur even with distinct non-zero numbers.
    F) False. The simple greedy algorithm is not optimal.
    G) False. A counterexample was found.
    H) True. An example is a set like {3,-3,5,-5} whose optimal solution is alternating.
    I) False. Minimizing max deviation is not the same as minimizing variance.
    J) True. The structure of optimal solutions for such sets tends to keep pairs adjacent to cancel out their sums, with at most one pair split.
    """
    true_statements = ['C', 'H', 'J']
    print("['C', 'H', 'J']")

solve()