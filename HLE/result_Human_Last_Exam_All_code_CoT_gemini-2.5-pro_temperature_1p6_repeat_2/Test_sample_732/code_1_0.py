import collections

def check_c4_free(num_vertices, edges):
    """
    Checks if a graph is C4-free.

    A graph is C4-free if and only if no two distinct vertices
    have more than one common neighbor.

    Args:
        num_vertices (int): The number of vertices in the graph.
        edges (list of tuples): The list of edges in the graph.

    Returns:
        bool: True if the graph is C4-free, False otherwise.
    """
    if num_vertices < 4:
        return True

    adj = collections.defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    for i in range(num_vertices):
        for j in range(i + 1, num_vertices):
            # Find common neighbors for vertices i and j
            # Using sets for efficient intersection
            neighbors_i = set(adj[i])
            neighbors_j = set(adj[j])
            
            common_neighbors = neighbors_i.intersection(neighbors_j)
            
            if len(common_neighbors) > 1:
                # Found two vertices with more than one common neighbor, so a C4 exists
                return False
                
    return True

def main():
    """
    Main function to construct the graph and find the maximum number of edges.
    """
    num_vertices = 8
    
    # Construct a C4-free graph with 8 vertices.
    # Start with a star graph K_1,7 (center 0, leaves 1-7)
    # Then add a matching on the leaves.
    edges = [
        (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), # Star graph edges
        (1, 2), (3, 4), (5, 6)                                  # Matching edges
    ]
    
    number_of_edges = len(edges)

    is_free = check_c4_free(num_vertices, edges)
    
    if is_free:
        print(f"A C4-free graph with {num_vertices} vertices has been constructed.")
        print("The edges are:")
        print(edges)
        
        # Calculate degrees to show the final equation
        degrees = collections.defaultdict(int)
        for u, v in edges:
            degrees[u] += 1
            degrees[v] += 1
        
        sum_of_degrees = sum(degrees.values())
        
        print("\nThe sum of the degrees of all vertices is:")
        degree_str = " + ".join(map(str, degrees.values()))
        print(f"{degree_str} = {sum_of_degrees}")
        
        print("\nThe number of edges is half the sum of the degrees:")
        print(f"Edges = {sum_of_degrees} / 2 = {number_of_edges}")
        
        print("\nBased on known results in extremal graph theory (Turán's problem for C4),")
        print(f"the maximum number of edges in a C4-free graph on {num_vertices} vertices is {number_of_edges}.")
        
    else:
        print("The constructed graph is not C4-free.")


if __name__ == "__main__":
    main()
