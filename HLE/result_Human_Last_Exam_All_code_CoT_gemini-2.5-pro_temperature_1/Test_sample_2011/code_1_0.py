import math

def solve_max_clique_sizes():
    """
    Calculates the maximum possible number of different clique sizes that can
    simultaneously appear as disjoint induced subgraphs in a graph on n vertices.
    This is equivalent to finding the maximum number of maximal cliques of
    distinct sizes.
    """
    n = 128

    print(f"The problem is to find the maximum number of different clique sizes in a graph with n = {n} vertices.")
    print("Under the interpretation that we are creating disjoint cliques of distinct sizes, we need to find the largest 'm' such that the sum of the first 'm' integers is at most n.")
    print("\nThe equation to solve is:")
    print(f"1 + 2 + ... + m <= {n}")
    print("Using the formula for the sum of the first 'm' integers, we get:")
    print(f"m * (m + 1) / 2 <= {n}")
    
    # We can multiply by 2 to simplify the inequality for integer calculations.
    rhs = 2 * n
    print(f"m * (m + 1) <= {rhs}")
    print("\nWe need to find the largest integer m that satisfies this inequality.")

    # We can solve for m by finding the root of m^2 + m - rhs = 0
    # m = (-1 + sqrt(1 + 4*rhs)) / 2
    m_float = (-1 + math.sqrt(1 + 4 * rhs)) / 2
    max_m = math.floor(m_float)
    
    print("\nLet's verify the result:")
    m_test_ok = max_m
    val_ok = m_test_ok * (m_test_ok + 1)
    print(f"For m = {m_test_ok}: {m_test_ok} * ({m_test_ok} + 1) = {val_ok}. Since {val_ok} <= {rhs}, this is a valid number of clique sizes.")
    
    m_test_fail = max_m + 1
    val_fail = m_test_fail * (m_test_fail + 1)
    print(f"For m = {m_test_fail}: {m_test_fail} * ({m_test_fail} + 1) = {val_fail}. Since {val_fail} > {rhs}, this is not possible.")

    print(f"\nThe maximum possible number of different clique sizes is {max_m}.")

solve_max_clique_sizes()