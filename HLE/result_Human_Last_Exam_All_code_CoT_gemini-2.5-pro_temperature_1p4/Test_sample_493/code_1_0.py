import numpy as np
from scipy.spatial import cKDTree

def calculate_average_constellation_size():
    """
    Simulates the constellation model to find the average number of stars per constellation.
    """
    # Number of stars to simulate. A larger number gives a more accurate estimate.
    N = 50000

    # 1. Generate N random star positions in a 2D square.
    # The specific size of the box does not affect the result due to uniform distribution.
    box_size = 1.0
    try:
        # Use a fixed seed for reproducibility of the output.
        rng = np.random.default_rng(seed=42)
        points = rng.random((N, 2)) * box_size
    except TypeError:
        # Fallback for older numpy versions that don't have default_rng
        np.random.seed(42)
        points = np.random.rand(N, 2) * box_size


    # 2. Find the nearest neighbor for each star.
    # We use a k-d tree for efficient nearest-neighbor searching.
    # The 'boxsize' parameter handles toroidal (wrap-around) boundaries to avoid edge effects.
    tree = cKDTree(points, boxsize=[box_size, box_size])
    # We query for the 2 nearest neighbors because the closest point to any point is itself.
    # The second-closest is its actual nearest neighbor.
    distances, indices = tree.query(points, k=2)
    # Get the indices of the nearest neighbors for each point.
    nn_indices = indices[:, 1]

    # 3. Build the graph representing the constellations.
    # An edge exists between star A and B if B is A's nearest neighbor (or vice versa).
    # We use an adjacency list to represent this undirected graph.
    adj = [[] for _ in range(N)]
    for i in range(N):
        j = nn_indices[i]
        # Add an undirected edge between star i and its nearest neighbor j
        adj[i].append(j)
        adj[j].append(i)

    # 4. Count the number of connected components (constellations).
    # We use a Breadth-First Search (BFS) to find all stars in each component.
    visited = [False] * N
    num_components = 0
    for i in range(N):
        if not visited[i]:
            # Found a new, unvisited star, so it must belong to a new component.
            num_components += 1
            q = [i]  # Start a BFS from this star to find all connected stars.
            visited[i] = True
            head = 0
            while head < len(q):
                u = q[head]
                head += 1
                # Add all unvisited neighbors to the queue to visit them later.
                for v in adj[u]:
                    if not visited[v]:
                        visited[v] = True
                        q.append(v)

    # 5. Calculate and print the average constellation size.
    # The average size is the total number of stars divided by the number of components.
    if num_components > 0:
        average_size = N / num_components
        print("Based on the simulation:")
        print(f"Total number of stars = {N}")
        print(f"Number of constellations (connected components) = {num_components}")
        print("\nAverage number of stars per constellation:")
        # The final equation with the numbers from the simulation
        print(f"{N} / {num_components} = {average_size:.4f}")
    else:
        print("No components found in the simulation.")

calculate_average_constellation_size()
<<<3.4228>>>