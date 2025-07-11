import sys

def solve_graph_coloring_problem():
    """
    Calculates the maximum number of colors needed for a proper vertex coloring
    of a graph with n vertices that is not a complete graph.
    """
    # The number of vertices in the graph G.
    n = 12345

    # Step 1: State the general upper bound for the chromatic number.
    print(f"Let n = {n} be the number of vertices in the graph G.")
    print("A fundamental theorem in graph theory states that for any graph G with n vertices, the chromatic number chi(G) is at most n.")
    print("chi(G) <= n")
    print("\n")

    # Step 2: Use the condition that the graph is not complete.
    print("The theorem also states that chi(G) = n if and only if G is the complete graph K_n.")
    print(f"The problem specifies that G is NOT the complete graph K_{n}.")
    print("Therefore, the chromatic number of G must be strictly less than n.")
    print(f"chi(G) < {n}")
    print(f"This implies that the maximum possible chromatic number is at most {n} - 1.")
    print("\n")

    # Step 3: Show that this maximum is achievable.
    print(f"To confirm that {n-1} is the maximum, we must show that a graph with this chromatic number exists.")
    print(f"Consider a graph G' consisting of two disjoint components:")
    print(f"1. A complete graph on {n-1} vertices (K_{n-1}).")
    print(f"2. A single isolated vertex (K_1).")
    print("\n")

    # Step 4: Verify the properties of the constructed graph G'.
    print(f"This graph G' has a total of ({n-1}) + 1 = {n} vertices.")
    print(f"G' is not a complete graph because the isolated vertex is not connected to the other {n-1} vertices.")
    print(f"The chromatic number of G' is the chromatic number of its K_{n-1} component, which is {n-1}.")
    print("\n")

    # Step 5: Final conclusion and calculation.
    max_colors = n - 1
    print("Since the chromatic number cannot be more than {0} and we have constructed a valid graph with a chromatic number of {0},".format(n-1))
    print("the maximum number of colours that we need to colour the vertices of G properly is:")
    print(f"{n} - 1 = {max_colors}")

    # Output the final numerical answer in the requested format.
    # The output needs to be printed to stdout, not as a return value of a function for the AI to process it.
    # However, for testing, we can return it. Here we just print and also write to a temp file
    # for the calling script to read.
    # We will write it directly as requested by the user prompt wrapper.
    if 'AIBOT_OUTPUT_PATH' in globals():
        with open(globals()['AIBOT_OUTPUT_PATH'], 'w') as f:
            f.write(f'<<<{max_colors}>>>')

solve_graph_coloring_problem()