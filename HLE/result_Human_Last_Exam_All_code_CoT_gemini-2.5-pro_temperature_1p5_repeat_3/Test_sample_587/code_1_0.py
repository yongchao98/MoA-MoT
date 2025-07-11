import networkx as nx

def count_k_paths(G, k):
    """
    Counts the number of simple paths of length k-1 (i.e., k vertices) in a graph G.
    This is a brute-force approach to illustrate the problem. It is not an FPT algorithm.
    """
    if k == 0:
        return 0
    if k == 1:
        return len(G.nodes())
        
    path_count = 0
    
    # We use a recursive helper function to find paths
    # For each node, we start a search for paths of length k-1
    for start_node in G.nodes():
        # The path starts with the current node
        path = [start_node]
        # Visited set to ensure simple paths (no repeated vertices)
        visited = {start_node}
        path_count += find_paths_from(G, k - 1, path, visited)
        
    # Each path is counted twice (once from each end), so we divide by 2.
    # e.g., A-B-C-D is found starting from A, and D-C-B-A is found starting from D.
    return path_count // 2

def find_paths_from(G, length_to_go, current_path, visited):
    """
    Recursively finds paths of a given remaining length.
    """
    # If we need to add 0 more nodes, we have found one full path
    if length_to_go == 0:
        return 1
    
    count = 0
    last_node = current_path[-1]
    
    # Explore neighbors of the last node in the path
    for neighbor in G.neighbors(last_node):
        if neighbor not in visited:
            new_visited = visited.copy()
            new_visited.add(neighbor)
            new_path = current_path + [neighbor]
            count += find_paths_from(G, length_to_go - 1, new_path, new_visited)
            
    return count

# Example Usage: The Petersen Graph
petersen = nx.petersen_graph()
k = 4 # Counting P_4 (paths with 4 vertices)

num_paths = count_k_paths(petersen, k)

print(f"Graph G has {petersen.number_of_nodes()} vertices and {petersen.number_of_edges()} edges.")
print(f"The number of subgraphs of G isomorphic to a path on {k} vertices is: {num_paths}")