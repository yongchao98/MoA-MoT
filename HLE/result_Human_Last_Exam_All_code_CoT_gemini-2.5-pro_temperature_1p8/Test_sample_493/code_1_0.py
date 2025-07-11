import numpy as np
from scipy.spatial import cKDTree
from collections import deque

def solve_constellation_problem():
    """
    This function simulates the star constellation model to find the average
    number of stars per constellation.

    A constellation is a connected component in a graph where each star is
    a node and is connected to its nearest neighbor.
    """
    # Number of stars to simulate. A larger number provides a more accurate result.
    N = 10000

    # 1. Generate N stars as random points in a 2D unit square.
    points = np.random.rand(N, 2)

    # 2. For each star, find its nearest neighbor.
    # We use scipy's cKDTree for an efficient search (O(N log N)).
    # We query for k=2 because the closest point to any point is itself.
    kdtree = cKDTree(points)
    # The second column [:, 1] gives the index of the nearest neighbor.
    _, nearest_neighbor_indices = kdtree.query(points, k=2)
    nearest_neighbors = nearest_neighbor_indices[:, 1]

    # 3. Build the graph.
    # We use an adjacency list to represent the connections. For finding
    # connected components, the graph is treated as undirected.
    adj_list = [[] for _ in range(N)]
    for i in range(N):
        neighbor = nearest_neighbors[i]
        # An edge exists if one star is the nearest neighbor of the other.
        adj_list[i].append(neighbor)
        adj_list[neighbor].append(i)

    # 4. Find all connected components (constellations) using BFS.
    visited = [False] * N
    num_components = 0
    
    for i in range(N):
        if not visited[i]:
            num_components += 1
            # Start a BFS traversal for the new component.
            q = deque([i])
            visited[i] = True
            while q:
                u = q.popleft()
                for v in adj_list[u]:
                    if not visited[v]:
                        visited[v] = True
                        q.append(v)

    # 5. Calculate and print the average constellation size.
    # The average size is the total number of stars divided by the number of constellations.
    average_size = N / num_components
    
    print(f"Simulation Parameters:")
    print(f"Total number of stars (N) = {N}")
    print("\nResults:")
    print(f"Number of constellations found = {num_components}")
    print("\nFinal Calculation:")
    print(f"Average stars per constellation = {N} / {num_components} = {average_size}")


if __name__ == "__main__":
    solve_constellation_problem()
