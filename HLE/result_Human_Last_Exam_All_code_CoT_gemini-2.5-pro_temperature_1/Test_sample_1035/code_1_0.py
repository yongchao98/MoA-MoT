import math

def solve_graph_problem():
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering all edges of the complete graph K_n for n=35.
    """
    n = 35

    print(f"The problem is to find the minimum total number of vertices in a family of bipartite graphs covering all edges of K_n for n = {n}.")
    print("This value is given by the formula: S = sum_{i=1 to n} ceil(log2(i)).")
    print("\nWe calculate this sum by grouping terms with the same value of ceil(log2(i)):")
    
    equation_parts = []
    total_sum = 0
    
    if n <= 0:
        print("n must be a positive integer.")
        return

    # Using a more robust way to group terms
    last_i = 0
    k_max = math.ceil(math.log2(n)) if n > 0 else 0

    for k in range(k_max + 1):
        # The range of i for which ceil(log2(i)) == k is (2^(k-1), 2^k].
        # For k=0, this is not well-defined, but we know ceil(log2(i))=0 only for i=1.
        if k == 0:
            if n >= 1:
                num_terms = 1
                contribution = num_terms * k
                total_sum += contribution
                equation_parts.append(f"{num_terms} * {k}")
            continue

        # For k >= 1
        lower_bound = 2**(k - 1)
        upper_bound = 2**k
        
        # Number of integers i in (lower_bound, upper_bound] that are also <= n
        count = max(0, min(n, upper_bound) - lower_bound)

        if count > 0:
            contribution = count * k
            total_sum += contribution
            equation_parts.append(f"{count} * {k}")

    final_equation = " + ".join(equation_parts)
    print(f"\nThe final calculation is: {final_equation} = {total_sum}")

solve_graph_problem()