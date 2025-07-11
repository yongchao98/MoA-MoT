import math

def solve_bipartite_covering():
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering all the edges of the complete graph K_n.
    """
    n = 35

    # The formula for the minimum number of vertices is n * ceil(log2(n)).
    
    # Step 1: Explain the value of n.
    print(f"The number of vertices in the complete graph is n = {n}.")
    
    # Step 2: Calculate ceil(log2(n)).
    # We need to find the smallest integer k such that 2^k >= n.
    # 2^5 = 32
    # 2^6 = 64
    # Since 32 < 35 <= 64, the ceiling of log2(35) is 6.
    log_val = math.log2(n)
    k = math.ceil(log_val)
    
    print(f"To find the minimum number of bipartite graphs for the construction, we calculate ceil(log2(n)).")
    print(f"log2({n}) is approximately {log_val:.4f}.")
    print(f"The ceiling of {log_val:.4f} is {k}.")
    
    # Step 3: Calculate the final result.
    result = n * k
    
    # Step 4: Print the final equation and result.
    print("\nThe minimum total number of vertices is given by the formula n * ceil(log2(n)).")
    print("The final calculation is:")
    print(f"{n} * {k} = {result}")

solve_bipartite_covering()