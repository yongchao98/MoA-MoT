import sys

def solve_clique_problem():
    """
    Calculates the maximum number of different clique sizes that can
    simultaneously appear as vertex-disjoint induced subgraphs in a
    graph of n vertices.

    This is equivalent to finding the largest integer m such that the sum of
    the first m positive integers is less than or equal to n.
    """
    # Set the total number of vertices for the graph.
    n = 128

    clique_sizes = []
    total_vertices_used = 0
    current_size = 1

    # Loop to find the largest number of cliques (m) with distinct sizes (1, 2, ..., m)
    # that can be formed without exceeding the total number of vertices.
    while total_vertices_used + current_size <= n:
        total_vertices_used += current_size
        clique_sizes.append(current_size)
        current_size += 1

    # The result is the number of different clique sizes we could fit.
    max_num_sizes = len(clique_sizes)

    # Build the equation string to display the calculation.
    # The numbers in the equation are the elements of the clique_sizes list.
    equation_str = " + ".join(map(str, clique_sizes))

    # Print the explanation and the final result including the equation.
    print(f"For a graph with n = {n} vertices:")
    print(f"We find the maximum number of disjoint cliques with distinct sizes (1, 2, 3, ...).")
    
    print("\nThe clique sizes that can be formed are:")
    # Using sys.stdout.write to prevent python interpreter from possibly wrapping long lines
    sys.stdout.write(", ".join(map(str, clique_sizes)) + "\n")

    print("\nThe final equation for the total vertices used is:")
    # Output each number in the final equation as requested.
    print(f"{equation_str} = {total_vertices_used}")
    
    print(f"\nThis construction uses {total_vertices_used} vertices, leaving {n - total_vertices_used} unused.")
    
    print(f"\nBased on this construction, the maximum possible number of different clique sizes is {max_num_sizes}.")

solve_clique_problem()