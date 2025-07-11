import math

def solve_graph_problem():
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering all the edges of the complete graph K_n.
    """
    n = 35

    # The problem asks for the minimum sum of the number of vertices in a family
    # of bipartite graphs that cover all edges of K_n.
    # A known theorem states that this minimum is given by the formula: n * ceil(log2(n)).

    # Step 1: We have n = 35.

    # Step 2: Calculate ceil(log2(n)).
    # We need to find the smallest integer k such that 2^k >= n.
    # 2^5 = 32
    # 2^6 = 64
    # Since 32 < 35 <= 64, the ceiling of log2(35) is 6.
    k = math.ceil(math.log2(n))

    # Step 3: Multiply n by the result from Step 2.
    min_total_vertices = n * k

    # Output the steps of the calculation and the final result.
    print(f"The problem is to find the minimum total vertices for K_n where n = {n}.")
    print(f"The formula for the minimum is n * ceil(log2(n)).")
    print(f"First, we calculate ceil(log2({n})).")
    print(f"log2({n}) is approximately {math.log2(n):.4f}.")
    print(f"The ceiling of {math.log2(n):.4f} is {k}.")
    print(f"Finally, we compute the result of the equation:")
    print(f"{n} * {k} = {min_total_vertices}")
    print(f"The minimum number of vertices is {min_total_vertices}.")

solve_graph_problem()