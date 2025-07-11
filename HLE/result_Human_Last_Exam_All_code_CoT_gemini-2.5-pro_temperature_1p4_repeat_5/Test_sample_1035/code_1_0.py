import math

def solve_bipartite_covering():
    """
    Calculates the minimum number of vertices in a family of bipartite graphs
    covering all edges of the complete graph K_n.
    """
    n = 35

    # The formula for the minimum total number of vertices is n * ceil(log2(n)).
    
    # Step 1: Explain the problem and formula
    print(f"The problem is to find the minimum total number of vertices in a family of bipartite graphs that cover all the edges of the complete graph K_n for n = {n}.")
    print("The formula for this minimum value is: n * ceil(log2(n))")
    print("-" * 20)
    
    # Step 2: Calculate log2(n)
    log2_n = math.log2(n)
    print(f"First, we calculate log2(n):")
    print(f"log2({n}) ≈ {log2_n:.4f}")
    
    # Step 3: Calculate ceil(log2(n))
    ceil_log2_n = math.ceil(log2_n)
    print(f"\nNext, we find the ceiling of this value:")
    print(f"ceil({log2_n:.4f}) = {ceil_log2_n}")
    
    # Step 4: Compute the final result
    min_vertices = n * ceil_log2_n
    print("\nFinally, we compute the total minimum number of vertices using the formula:")
    # As requested, output each number in the final equation.
    print(f"{n} * {ceil_log2_n} = {min_vertices}")

if __name__ == '__main__':
    solve_bipartite_covering()