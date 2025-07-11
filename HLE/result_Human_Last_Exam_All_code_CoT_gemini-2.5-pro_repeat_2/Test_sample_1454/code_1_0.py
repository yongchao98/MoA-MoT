import collections

def solve():
    """
    Calculates the number of connected components for the fractal F.

    The number of components is determined by the connectivity of a graph
    representing the iterated function system.
    
    The vertices of the graph are the vectors in D. An edge exists between
    d1 = (x1, y1) and d2 = (x2, y2) if the generated squares overlap, which
    occurs if |x1 - x2| <= 1 and |y1 - y2| <= 1.
    """
    D = [(0,0), (0,1), (0,2), (0,3), (3,0), (3,1), (3,2), (3,3)]
    
    # Adjacency list representation of the graph
    adj = collections.defaultdict(list)
    
    # Build the graph
    for i in range(len(D)):
        for j in range(i + 1, len(D)):
            d1 = D[i]
            d2 = D[j]
            
            # Check for edge condition
            if abs(d1[0] - d2[0]) <= 1 and abs(d1[1] - d2[1]) <= 1:
                adj[i].append(j)
                adj[j].append(i)

    # Count connected components using graph traversal (BFS)
    visited = set()
    num_components = 0
    
    for i in range(len(D)):
        if i not in visited:
            num_components += 1
            q = collections.deque([i])
            visited.add(i)
            component_nodes = []
            while q:
                u = q.popleft()
                component_nodes.append(D[u])
                for v in adj[u]:
                    if v not in visited:
                        visited.add(v)
                        q.append(v)
            # The final equation requires outputting each number.
            # We interpret this as showing the components found.
            print(f"Component {num_components}: {component_nodes}")

    print(f"\nThe smallest possible number of nondegenerate and locally connected components of F is equal to the number of connected components of the IFS graph.")
    print(f"Final equation: Number of components = {num_components}")


solve()
