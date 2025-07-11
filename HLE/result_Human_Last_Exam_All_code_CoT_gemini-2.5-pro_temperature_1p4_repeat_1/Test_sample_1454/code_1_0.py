import collections

def solve():
    """
    Solves for the number of components by analyzing the connectivity graph of the IFS.
    """
    
    # The set D defines the translations for the IFS.
    D = [(0,0), (0,1), (0,2),(0,3),(3,0), (3,1), (3,2),(3,3)]
    
    # Generate the 8 squares (represented by their intervals)
    # S_i = [x_min, x_max, y_min, y_max]
    squares = []
    for d in D:
        d1, d2 = d
        squares.append([d1 / 4.0, (d1 + 1) / 4.0, d2 / 4.0, (d2 + 1) / 4.0])

    # Build adjacency list for the graph
    adj = collections.defaultdict(list)
    num_squares = len(squares)

    for i in range(num_squares):
        for j in range(i + 1, num_squares):
            # Check for intersection between square i and square j
            s1 = squares[i]
            s2 = squares[j]
            
            # x_overlap is true if [s1[0], s1[1]] and [s2[0], s2[1]] intersect
            x_overlap = (s1[0] <= s2[1] and s1[1] >= s2[0])
            
            # y_overlap is true if [s1[2], s1[3]] and [s2[2], s2[3]] intersect
            y_overlap = (s1[2] <= s2[3] and s1[3] >= s2[2])
            
            if x_overlap and y_overlap:
                adj[i].append(j)
                adj[j].append(i)

    # Find number of connected components using Breadth-First Search (BFS)
    visited = [False] * num_squares
    num_components = 0
    for i in range(num_squares):
        if not visited[i]:
            num_components += 1
            q = collections.deque([i])
            visited[i] = True
            while q:
                u = q.popleft()
                for v in adj[u]:
                    if not visited[v]:
                        visited[v] = True
                        q.append(v)
    
    print(f"The number of components in the contact graph is: {num_components}")

solve()