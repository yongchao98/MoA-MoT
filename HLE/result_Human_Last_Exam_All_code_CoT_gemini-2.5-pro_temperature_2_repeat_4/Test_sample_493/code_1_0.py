import numpy as np
from scipy.spatial import cKDTree
from collections import deque

def calculate_average_constellation_size():
    """
    Simulates the constellation model and calculates the average number of stars
    per constellation.
    """
    # Step 1: Set up simulation parameters. A large number of stars gives a
    # more accurate average, converging to the theoretical value.
    N = 50000  # Total number of stars

    # Step 2: Generate N random stars in a 2D unit square.
    stars = np.random.rand(N, 2)

    # Step 3: Find the nearest neighbor for each star.
    # We use a k-d tree for efficient O(N log N) searching.
    # We query for k=2 because the closest point to any given star is the star
    # itself. Its nearest neighbor will be the second-closest point.
    tree = cKDTree(stars)
    # The 'indices' array will store the index of the neighbors.
    # indices[:, 0] will be the point itself, indices[:, 1] will be the nearest neighbor.
    distances, indices = tree.query(stars, k=2, workers=-1)
    
    # Extract the nearest neighbor indices for each star.
    nearest_neighbor_indices = indices[:, 1]

    # Step 4: Build the undirected graph representing the constellations.
    # We use an adjacency list, where adj[i] stores the neighbors of star i.
    # We use sets initially to automatically handle duplicate edges.
    adj_sets = [set() for _ in range(N)]
    for i in range(N):
        j = nearest_neighbor_indices[i]
        # Add an undirected edge between star i and its nearest neighbor j.
        adj_sets[i].add(j)
        adj_sets[j].add(i)

    # Convert sets to lists for standard traversal algorithms.
    adj_list = [list(s) for s in adj_sets]

    # Step 5: Count the number of connected components (constellations) using BFS.
    visited = [False] * N
    num_components = 0
    for i in range(N):
        if not visited[i]:
            # This star is part of a new, unvisited constellation.
            num_components += 1
            # Start a Breadth-First Search (BFS) to find all stars in this component.
            q = deque([i])
            visited[i] = True
            while q:
                u = q.popleft()
                for v in adj_list[u]:
                    if not visited[v]:
                        visited[v] = True
                        q.append(v)

    # Step 6: Calculate the average size and print the result.
    if num_components > 0:
        average_size = N / num_components
    else:
        average_size = 0 # Should not happen for N > 0

    print("--- Simulation Results ---")
    print(f"Total stars simulated (N): {N}")
    print(f"Number of constellations found: {num_components}")
    print("\nCalculation for the average number of stars per constellation:")
    print(f"{N} / {num_components} = {average_size:.4f}")

if __name__ == '__main__':
    calculate_average_constellation_size()