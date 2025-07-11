import numpy as np
from collections import deque
# The scipy library is used for its efficient k-d tree implementation,
# which allows for fast nearest neighbor searching.
# If you don't have it, install it via: pip install scipy
from scipy.spatial import cKDTree

def calculate_average_constellation_size():
    """
    Simulates star constellations to find the average number of stars per constellation.
    """
    # Use a fixed seed for reproducibility of the random point generation.
    np.random.seed(42)

    # --- Step 1: Model the stars as N random points in a 2D plane ---
    # The average constellation size should be independent of N for large N.
    N = 10000  # Number of stars

    # Generate N random points in a 2D unit square.
    points = np.random.rand(N, 2)

    # --- Step 2: Form connections by finding the nearest neighbor for each star ---
    # Use a k-d tree for efficient nearest neighbor search (O(N log N)).
    # We query for k=2 because the closest point to any given point is itself.
    tree = cKDTree(points)
    # The `indices` array will contain the index of the nearest neighbor for each point.
    # indices[:, 0] will be the index of the point itself.
    # indices[:, 1] will be the index of its true nearest neighbor.
    _, nn_indices = tree.query(points, k=2)
    nearest_neighbors = nn_indices[:, 1]

    # Create an adjacency list for the undirected graph.
    # For each star `i`, an edge is added to its nearest neighbor `j`.
    adj = [[] for _ in range(N)]
    for i in range(N):
        j = nearest_neighbors[i]
        # Add an undirected edge {i, j} to the graph representation.
        adj[i].append(j)
        adj[j].append(i)

    # --- Step 3: Identify constellations (connected components) ---
    # We use a breadth-first search (BFS) to find all nodes connected to a
    # starting node, marking them as visited. We repeat this for all unvisited
    # nodes. The number of times we start a new traversal is the number of
    # connected components.
    visited = [False] * N
    num_components = 0
    for i in range(N):
        if not visited[i]:
            num_components += 1
            q = deque([i])
            visited[i] = True
            while q:
                u = q.popleft()
                for v in adj[u]:
                    if not visited[v]:
                        visited[v] = True
                        q.append(v)

    # --- Step 4: Calculate the average size of a constellation ---
    average_size = N / num_components

    # --- Output the results ---
    print("This simulation models stars as points uniformly distributed in a 2D plane.")
    print("Each star is connected by an edge to its single nearest neighbor, forming a graph.")
    print("The 'constellations' are the connected components of this graph.")
    print("\n--- Simulation Results ---")
    print(f"Total number of stars simulated: {N}")
    print(f"Number of separate constellations found: {num_components}")
    print("\n--- Final Calculation ---")
    print(f"Average Stars per Constellation = Total Stars / Number of Constellations")
    print(f"                                = {N} / {num_components} = {average_size:.4f}")

if __name__ == '__main__':
    calculate_average_constellation_size()