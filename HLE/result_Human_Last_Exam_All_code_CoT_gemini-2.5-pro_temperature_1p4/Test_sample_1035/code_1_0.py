import math

def solve_graph_problem():
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering all the edges of the complete graph K_n for n = 35.

    This problem is a known result in graph theory, solved by Fan Chung Graham.
    The formula for the minimum sum of vertices is n * ceil(log2(n)).
    """
    n = 35

    # Step 1: Calculate log base 2 of n.
    log_n_base_2 = math.log2(n)

    # Step 2: Find the ceiling of the result.
    # 2^5 = 32 and 2^6 = 64, so log2(35) is between 5 and 6.
    # The ceiling will be 6.
    ceil_log_n = math.ceil(log_n_base_2)

    # Step 3: Multiply by n to get the final answer.
    min_total_vertices = n * ceil_log_n

    # Output the explanation and the breakdown of the calculation.
    print(f"The problem is to find the minimum total number of vertices in a family of bipartite graphs covering K_n for n = {n}.")
    print("According to a theorem by Fan Chung Graham, this value is given by the formula: n * ceil(log2(n)).")
    print("\nCalculation steps:")
    print(f"1. n = {n}")
    print(f"2. log2({n}) is approximately {log_n_base_2:.4f}.")
    print(f"3. The ceiling of log2({n}) is {ceil_log_n}.")
    print("\nThe final equation with the numbers plugged in is:")
    print(f"{n} * {ceil_log_n} = {min_total_vertices}")

solve_graph_problem()