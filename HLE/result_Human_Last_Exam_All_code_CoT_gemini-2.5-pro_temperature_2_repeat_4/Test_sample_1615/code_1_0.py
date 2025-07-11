import sys

def solve_coloring_problem():
    """
    Calculates the maximum number of colors needed to properly color a graph G
    with 12345 vertices, given that G is not a complete graph.
    """
    
    # The number of vertices in the graph G
    n = 12345
    
    print(f"The graph G has n = {n} vertices.")
    print("The problem states that G is not a complete graph (K_n).")
    print("-" * 30)

    # Step 1: Upper bound for the chromatic number (χ)
    print("Step 1: Finding the upper bound for the number of colors.")
    print(f"For any graph with {n} vertices, the chromatic number χ is at most {n}.")
    print(f"A graph has a chromatic number χ = {n} if and only if it is the complete graph K_{n}.")
    print(f"Since G is not a complete graph, its chromatic number must be less than {n}.")
    print(f"Therefore, χ(G) ≤ {n} - 1.")
    print("-" * 30)
    
    # Step 2: Show that the maximum is achievable
    print("Step 2: Proving the maximum is achievable.")
    print("We need to show that a non-complete graph with a chromatic number of n - 1 can exist.")
    print(f"Let's construct such a graph: Consider a complete graph on {n-1} vertices, K_{n-1}, and add one single isolated vertex.")
    print(f"This new graph has ({n-1}) + 1 = {n} vertices.")
    print("It is not a complete graph because the isolated vertex is not connected to any other vertex.")
    print(f"The K_{n-1} component requires {n-1} colors for a proper coloring.")
    print(f"The chromatic number of this constructed graph is therefore {n-1}.")
    print("-" * 30)
    
    # Step 3: Conclusion and Final Calculation
    print("Step 3: Conclusion.")
    print("Since the number of colors cannot exceed n-1, and we have found a graph that requires n-1 colors,")
    print("the maximum number of colors is n-1.")
    
    max_colors = n - 1
    
    # Final output as an equation
    print("\nFinal Answer:")
    print(f"The maximum number of colours is {n} - 1 = {max_colors}")

solve_coloring_problem()

# Flushing the output to be safe in some environments.
sys.stdout.flush()