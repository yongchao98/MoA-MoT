import collections

def solve():
    """
    Solves the problem by building and analyzing the intersection graph of the IFS.
    """
    
    # The set D defines the translations for the Iterated Function System (IFS).
    D = [(0, 0), (0, 1), (0, 2), (0, 3), (3, 0), (3, 1), (3, 2), (3, 3)]
    
    # The scaling factor in the equation.
    scaling_factor = 4
    
    print("Numbers from the equation:")
    print(f"Scaling factor: {scaling_factor}")
    print("Points in D:")
    for point in D:
        print(point)
    print("-" * 20)
    
    num_maps = len(D)
    adj = collections.defaultdict(list)
    
    # Build the adjacency list for the intersection graph.
    # An edge exists between d1 and d2 if their corresponding squares intersect.
    # f_d([0,1]^2) = [d_x/4, (d_x+1)/4] x [d_y/4, (d_y+1)/4]
    # Intersection occurs iff x-intervals overlap and y-intervals overlap.
    # x-intervals [0, 1/4] and [3/4, 1] only overlap if they are the same. So d1_x must equal d2_x.
    # y-intervals [y/4, (y+1)/4] overlap if |d1_y - d2_y| <= 1.
    for i in range(num_maps):
        for j in range(i + 1, num_maps):
            d1 = D[i]
            d2 = D[j]
            
            # Check for x-overlap
            if d1[0] == d2[0]:
                # Check for y-overlap
                if abs(d1[1] - d2[1]) <= 1:
                    adj[i].append(j)
                    adj[j].append(i)

    # Count the number of connected components using BFS.
    visited = [False] * num_maps
    num_components = 0
    for i in range(num_maps):
        if not visited[i]:
            num_components += 1
            queue = collections.deque([i])
            visited[i] = True
            while queue:
                u = queue.popleft()
                for v in adj[u]:
                    if not visited[v]:
                        visited[v] = True
                        queue.append(v)
                        
    print(f"The smallest possible number of nondegenerate and locally connected components is interpreted as the number of connected components of the IFS intersection graph.")
    print(f"Number of components: {num_components}")

solve()