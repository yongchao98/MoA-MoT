import networkx as nx

def count_k_connected_graphs(vertices, connectivity_k):
    """
    Counts the number of k-vertex-connected simple nonisomorphic graphs.
    
    Args:
        vertices (int): The number of vertices in the graphs.
        connectivity_k (int): The minimum vertex connectivity required.
    """
    
    print(f"Searching for simple nonisomorphic graphs with {vertices} vertices and vertex connectivity >= {connectivity_k}.")
    
    count = 0
    # Generate all nonisomorphic graphs with the specified number of vertices.
    graph_generator = nx.nonisomorphic_graphs(vertices)
    
    # Iterate through each unique graph structure.
    for G in graph_generator:
        # A graph must have at least k+1 vertices to be k-connected.
        # This condition is met since vertices > connectivity_k.
        
        # Calculate the vertex connectivity of the graph.
        # A disconnected graph has connectivity 0.
        # A connected graph with a cut vertex has connectivity 1.
        # A 2-connected graph has connectivity >= 2.
        # The complete graph K_n has connectivity n-1.
        if nx.node_connectivity(G) >= connectivity_k:
            count += 1
            
    return count

# --- Main execution ---
# The problem asks for 2-vertex-connected graphs with 5 vertices.
num_vertices = 5
required_connectivity = 2

# Calculate the result using the function
result = count_k_connected_graphs(num_vertices, required_connectivity)

# Output the final answer in a sentence, showing the numbers involved.
print("\n--- Final Answer ---")
print(f"The number of simple nonisomorphic graphs with {num_vertices} vertices that are {required_connectivity}-vertex-connected is: {result}")