import math

def solve_clique_problem():
    """
    This script solves for the maximum number of different clique sizes in a graph
    with n=128 vertices, interpreted as finding the maximum number of maximal
    cliques of distinct sizes.
    """
    n = 128

    # We need to find the largest integer 'm' such that 1 + 2 + ... + m <= n.
    # This corresponds to partitioning n into the maximum number of distinct parts.
    m = 0
    current_sum = 0
    
    # We iteratively add the next integer until the sum exceeds n.
    while current_sum + (m + 1) <= n:
        m += 1
        current_sum += m

    # The final 'm' is the maximum number of distinct parts.
    # The prompt requires printing the final equation.
    equation_numbers = list(range(1, m + 1))
    equation_str = " + ".join(map(str, equation_numbers))

    print(f"The problem is to find the largest integer 'm' such that the sum of the first 'm' positive integers is at most {n}.")
    print(f"This can be expressed by the inequality: m * (m + 1) / 2 <= {n}.")
    print(f"\nBy solving this, we find the maximum value for 'm' is {m}.")
    
    print("\nThe sum of the sizes for the first m cliques is:")
    print(f"{equation_str} = {current_sum}")
    
    # We can construct a graph with these properties. For example, a disjoint union
    # of cliques of sizes 1, 2, ..., 14, and a final clique using the remaining
    # vertices: 15 + (128 - 120) = 23. The sizes {1, 2,..., 14, 23} are 15 distinct sizes.
    remaining_vertices = n - current_sum
    print(f"\nThis sum of {current_sum} uses up most of the {n} vertices, with {remaining_vertices} remaining.")
    print(f"If we tried m = {m + 1}, the sum would be {current_sum + m + 1}, which exceeds {n}.")
    print(f"\nTherefore, the maximum possible number of different clique sizes is {m}.")

solve_clique_problem()