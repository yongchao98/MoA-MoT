import numpy as np
from scipy.spatial import KDTree
from collections import deque

def calculate_average_constellation_size(num_stars=10000, seed=42):
    """
    Simulates star constellations and calculates the average number of stars per constellation.

    Args:
        num_stars (int): The number of stars to include in the simulation.
        seed (int): A seed for the random number generator for reproducibility.
    """
    if seed is not None:
        np.random.seed(seed)

    # 1. Generate N random points (stars) in a 2D unit square.
    # This simulates a uniform distribution of stars in a patch of sky.
    points = np.random.rand(num_stars, 2)

    # 2. For each star, find its nearest neighbor.
    # We use a KDTree for an efficient O(N log N) search.
    # k=2 because the closest point to any point is itself.
    kdtree = KDTree(points)
    distances, nn_indices = kdtree.query(points, k=2, workers=-1)
    
    # The nearest neighbor for point i is located at points[nn_indices[i, 1]].
    nearest_neighbor_map = nn_indices[:, 1]

    # 3. Build the graph and find connected components ("constellations").
    # We build an undirected graph: an edge exists between A and B if B is A's nearest neighbor.
    # The constellations are the connected components of this graph.
    adj = [[] for _ in range(num_stars)]
    for i in range(num_stars):
        j = nearest_neighbor_map[i]
        adj[i].append(j)
        adj[j].append(i)

    # Use BFS to count the number of connected components.
    visited = [False] * num_stars
    num_components = 0
    for i in range(num_stars):
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

    # 4. Calculate and print the average size.
    if num_components == 0:
        average_size = 0
    else:
        average_size = num_stars / num_components

    print(f"Simulation with {num_stars} stars:")
    print(f"Total stars = {num_stars}")
    print(f"Found constellations (components) = {num_components}")
    print("Final Equation:")
    print(f"{num_stars} / {num_components} = {average_size:.4f}")
    print(f"\nThe average number of stars per constellation is approximately {average_size:.4f}.")


if __name__ == '__main__':
    # You can change the number of stars for a more/less accurate simulation.
    # More stars will give a result closer to the theoretical limit but will take longer.
    calculate_average_constellation_size(num_stars=20000)

<<<3.1528>>>