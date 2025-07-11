import math

def solve():
    """
    This function calculates the maximum possible number of different clique sizes
    that can simultaneously appear as induced subgraphs of a single graph with n=128 vertices.

    The problem is interpreted as finding the maximum number of different *maximal* clique sizes.
    To maximize this number, m, we choose the smallest possible sizes: 1, 2, ..., m.

    We use a vertex-efficient "chained clique" construction:
    - A disjoint K_1 (maximal clique of size 1).
    - A chain of cliques K_2, K_3, ..., K_m, where K_i and K_{i+1} share a vertex.

    The number of vertices for this construction is N(m) = 1 + m + (m-1)*(m-2)/2.
    We find the largest integer m for which N(m) <= 128.
    """
    n = 128
    
    m = 1
    while True:
        # Calculate the number of vertices needed for m maximal cliques of sizes 1...m
        if m == 1:
            # For m=1, we just need a K_1.
            vertices_needed = 1
        else:
            # For m > 1, use the chained-clique formula.
            # sum_of_additions is for the chain K_m, ..., K_2
            chain_vertices = m + (m - 2) * (m - 1) // 2
            # Total vertices = K_1 + chain
            vertices_needed = 1 + chain_vertices
        
        if vertices_needed > n:
            # This m is too large, so the previous value was the maximum.
            m = m - 1
            break
        
        m += 1
        
    final_m = m
    
    # Recalculate values for the final m to display the equation.
    if final_m == 1:
        final_vertices_needed = 1
        sum_val = 0
    else:
        sum_val = (final_m - 1) * (final_m - 2) // 2
        final_vertices_needed = 1 + final_m + sum_val

    print(f"The number of vertices is n = {n}.")
    print("We aim to find the maximum number (m) of distinct maximal clique sizes.")
    print("To maximize m, we choose clique sizes 1, 2, ..., m.")
    print("Using a chained-clique construction, the number of vertices required is:")
    print("N(m) = 1 + m + (m-1)*(m-2)/2")
    print("")
    print(f"For m = {final_m}:")
    print(f"The equation is: 1 + {final_m} + ({final_m}-1)*({final_m}-2)/2 = 1 + {final_m} + {sum_val} = {final_vertices_needed}")
    print(f"Since {final_vertices_needed} <= {n}, it is possible to have {final_m} different sizes.")
    print("")

    # Check the next value of m to confirm it's not possible.
    next_m = final_m + 1
    next_sum_val = (next_m - 1) * (next_m - 2) // 2
    next_vertices_needed = 1 + next_m + next_sum_val
    print(f"For m = {next_m}:")
    print(f"The equation is: 1 + {next_m} + ({next_m}-1)*({next_m}-2)/2 = 1 + {next_m} + {next_sum_val} = {next_vertices_needed}")
    print(f"Since {next_vertices_needed} > {n}, it is not possible to have {next_m} different sizes with this construction.")
    print("")
    print(f"Thus, the maximum possible number of different clique sizes is {final_m}.")

solve()
<<<16>>>