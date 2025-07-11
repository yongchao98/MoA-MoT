import collections

def count_fractal_components():
    """
    This function determines the number of components of the fractal set F
    by modeling the problem as finding connected components in a graph.
    """
    
    # The set of translation vectors defining the Iterated Function System.
    D = [(0, 0), (0, 1), (0, 2), (0, 3), (3, 0), (3, 1), (3, 2), (3, 3)]
    
    # Create a mapping from vector to a unique integer index for easier graph processing.
    nodes = {d: i for i, d in enumerate(D)}
    num_nodes = len(D)
    
    # Build an adjacency list representation of the connectivity graph.
    # An edge exists between d_i and d_j if their corresponding squares intersect.
    # This happens if they are in the same column (d_i.x == d_j.x)
    # and their y-indices are adjacent or the same (|d_i.y - d_j.y| <= 1).
    adj = collections.defaultdict(list)
    for i in range(num_nodes):
        for j in range(i, num_nodes):
            d_i = D[i]
            d_j = D[j]
            
            # Check for the connectivity condition.
            if d_i[0] == d_j[0] and abs(d_i[1] - d_j[1]) <= 1:
                # Add an edge between the nodes in the graph.
                adj[nodes[d_i]].append(nodes[d_j])
                adj[nodes[d_j]].append(nodes[d_i])

    print("Step 1: The fractal F is composed of scaled and translated copies of itself.")
    print("Its connectivity is determined by how these copies connect.")
    print("This can be modeled as a graph where the nodes are the transformations.")
    
    print("\nStep 2: Building and analyzing the graph based on the set D.")
    print(f"The set of vectors is D = {D}")
    print("The graph has two distinct, non-connecting groups of nodes:")
    print("Group 1 (left bar): All nodes with x-coordinate 0.")
    print("Group 2 (right bar): All nodes with x-coordinate 3.")

    # Use a graph traversal algorithm (like BFS) to count connected components.
    visited = set()
    num_components = 0
    
    for i in range(num_nodes):
        if i not in visited:
            num_components += 1
            # Start a traversal for a new component.
            q = collections.deque([i])
            visited.add(i)
            component_nodes = []
            while q:
                u = q.popleft()
                component_nodes.append(D[u])
                for v_idx in adj[u]:
                    if v_idx not in visited:
                        visited.add(v_idx)
                        q.append(v_idx)
            print(f"-> Found Component {num_components} formed by vectors: {sorted(component_nodes)}")
            
    print("\nStep 3: Calculating the total number of components.")
    print("Each group of vectors forms a single connected component in the graph.")
    print("The total number of components of F is the sum of these.")
    print("Number of components from Group 1: 1")
    print("Number of components from Group 2: 1")
    print(f"Final equation: 1 + 1 = {num_components}")
    print("These components are nondegenerate and locally connected fractal arcs.")

count_fractal_components()
<<<2>>>