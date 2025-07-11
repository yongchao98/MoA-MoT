import networkx as nx

def solve_graph_counting():
    """
    This function calculates the number of 2-vertex-connected simple
    nonisomorphic graphs with 5 vertices by iterating through the networkx
    graph atlas.
    """
    try:
        # Step 1: Get all nonisomorphic graphs on 5 vertices from the atlas.
        all_graphs_from_atlas = nx.graph_atlas_g()
        graphs_n5 = [g for g in all_graphs_from_atlas if g.number_of_vertices() == 5]
        total_graphs = len(graphs_n5)

        # Step 2: Categorize each graph by its vertex connectivity.
        # We will count graphs with connectivity 0 (disconnected), 1, and >=2.
        connectivity_counts = {0: 0, 1: 0, '2+': 0}

        for g in graphs_n5:
            # A graph must have more than one vertex to compute connectivity.
            # All our graphs have 5 vertices, so this is fine.
            if not nx.is_connected(g):
                # Disconnected graphs have vertex connectivity 0.
                connectivity_counts[0] += 1
            else:
                # For connected graphs, find the node connectivity.
                conn = nx.node_connectivity(g)
                if conn == 1:
                    # 1-vertex-connected (has a cut vertex).
                    connectivity_counts[1] += 1
                else:
                    # 2-vertex-connected or more (biconnected).
                    # For n=5, max connectivity is n-1=4.
                    connectivity_counts['2+'] += 1
        
        disconnected_count = connectivity_counts[0]
        one_connected_count = connectivity_counts[1]
        two_or_more_connected_count = connectivity_counts['2+']

        # Step 3: Print the results, including the "equation".
        print(f"Analysis of simple nonisomorphic graphs with 5 vertices:")
        print(f"Total number of graphs = {total_graphs}")
        print(f" - Number of disconnected graphs (connectivity 0) = {disconnected_count}")
        print(f" - Number of 1-connected graphs (separable) = {one_connected_count}")
        print(f" - Number of 2-vertex-connected graphs (biconnected) = {two_or_more_connected_count}")
        
        print("\nThe final count can be represented by the equation:")
        print(f"{total_graphs} (Total) - {disconnected_count} (Disconnected) - {one_connected_count} (1-Connected) = {two_or_more_connected_count}")
        
        print(f"\nThus, the number of 2-vertex-connected simple nonisomorphic graphs with 5 vertices is {two_or_more_connected_count}.")

    except ImportError:
        print("NetworkX library not found. Please install it using 'pip install networkx'")
    except Exception as e:
        print(f"An error occurred: {e}")

solve_graph_counting()