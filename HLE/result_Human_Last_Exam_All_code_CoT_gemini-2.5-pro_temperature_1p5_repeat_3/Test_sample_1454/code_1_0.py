def solve():
    """
    Calculates the number of connected components of the IFS graph.
    This is interpreted as the number of 'macro-components' of the fractal.
    """
    D = [(0,0), (0,1), (0,2),(0,3),(3,0), (3,1), (3,2),(3,3)]
    
    # Adjacency list representation of the graph
    adj = {d: [] for d in D}
    
    # Build the graph
    for i in range(len(D)):
        for j in range(i + 1, len(D)):
            d1 = D[i]
            d2 = D[j]
            
            # Check for edge condition
            # x-coordinates must be the same
            cond_x = (d1[0] == d2[0])
            # y-coordinates must be adjacent
            cond_y = (abs(d1[1] - d2[1]) <= 1)
            
            if cond_x and cond_y:
                adj[d1].append(d2)
                adj[d2].append(d1)

    # Count connected components using BFS/DFS
    num_components = 0
    visited = set()
    
    for d_start in D:
        if d_start not in visited:
            num_components += 1
            # Start a traversal (e.g., BFS) from this node
            q = [d_start]
            visited.add(d_start)
            while q:
                u = q.pop(0)
                for v in adj[u]:
                    if v not in visited:
                        visited.add(v)
                        q.append(v)
                        
    print(f"The set of vectors is D = {D}")
    print("A connection graph is built with vertices from D.")
    print("An edge connects d1 = (x1, y1) and d2 = (x2, y2) if their generated regions overlap.")
    print("This occurs if x1 = x2 and |y1 - y2| <= 1.")
    print("Let's find the number of connected components in this graph.")
    print(f"The graph has {len(D)} vertices and connections based on the rules above.")
    print(f"Number of connected components found: {num_components}")
    print("\nBased on this analysis, the number of components is:")
    # The final equation is trivial, but the instructions require it.
    final_result = num_components
    print(f"{final_result}")

solve()
<<<2>>>