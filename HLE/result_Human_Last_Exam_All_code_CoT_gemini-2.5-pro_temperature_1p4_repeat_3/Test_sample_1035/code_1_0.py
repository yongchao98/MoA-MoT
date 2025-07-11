import math

def calculate_min_vertices(n):
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering all edges of the complete graph K_n.

    This function iterates through possible numbers of bipartite graphs (k)
    and calculates the minimum cost for each k. The cost is the total number
    of non-'*' symbols in a separating system of n codes.
    """

    print(f"Finding the minimum number of vertices for K_{n} where n = {n}.\n")
    print("This is equivalent to a coding theory problem:")
    print("1. Assign a unique code vector from {0, 1, *}^k to each of the n vertices.")
    print("2. For any pair of vertices, their codes must be 'separated' (one has 0, the other 1 in some position).")
    print("3. The cost is the total number of non-'*' symbols. We want to minimize this cost.")
    print("4. We can vary k (the number of bipartite graphs) to find the global minimum cost.\n")

    min_total_vertices = float('inf')
    optimal_k = -1

    # We must have enough codes for n vertices. 3^k - 1 >= n.
    # For n=35, 3^3-1=26 < 35, 3^4-1=80 > 35. So k must be at least 4.
    # We check a reasonable range of k values.
    # The cost function C(k) = 2n-2k for k < n/2 and C(k)=n for k >= n/2 suggests the minimum is around n/2.
    k_start = math.ceil(math.log(n + 1, 3))
    k_end = n # A safe upper bound for k to check

    print(f"Calculating costs for different numbers of bipartite graphs (k) from {k_start} to {k_end}:\n")

    for k in range(k_start, k_end + 1):
        vertices_to_code = n
        current_total_vertices = 0
        
        print(f"For k = {k}:")
        
        calculation_explanation = []

        for j in range(1, k + 1): # j is the weight of the code
            if vertices_to_code == 0:
                break
            
            try:
                # Number of codes with weight j
                num_codes_at_weight_j = math.comb(k, j) * (2 ** j)
            except ValueError:
                num_codes_at_weight_j = 0
            
            if num_codes_at_weight_j > 0:
                codes_to_take = min(vertices_to_code, num_codes_at_weight_j)
                cost_at_this_weight = codes_to_take * j
                current_total_vertices += cost_at_this_weight
                vertices_to_code -= codes_to_take

                # Add explanation for this step
                num_needed = min(vertices_to_code + codes_to_take, n)
                explanation_line = f"  - Weight {j}: Taking {codes_to_take} codes. Cost added: {codes_to_take} * {j} = {cost_at_this_weight}"
                calculation_explanation.append(explanation_line)


        for line in calculation_explanation:
            print(line)
        
        print(f"  Total cost for k={k}: {current_total_vertices}\n")

        if current_total_vertices < min_total_vertices:
            min_total_vertices = current_total_vertices
            optimal_k = k

    print("--- Summary ---")
    print(f"The minimum total number of vertices is {min_total_vertices}.")
    print(f"This minimum is achieved with k = {optimal_k} bipartite graphs (and any k > {optimal_k}).")


if __name__ == '__main__':
    n = 35
    calculate_min_vertices(n)
    # The calculation shows the minimum is 35.
    print("The final answer is the minimum value found.")
    