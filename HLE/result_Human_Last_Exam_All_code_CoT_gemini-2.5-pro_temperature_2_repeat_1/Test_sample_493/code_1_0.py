import numpy as np
from scipy.spatial import cKDTree
from collections import deque

def solve_constellation_problem():
    """
    Simulates the constellation problem to find the average number of stars
    per constellation.
    """
    # The number of stars to use in the simulation.
    # A larger number provides a more accurate approximation of the theoretical value.
    num_stars = 20000

    print(f"Running simulation with {num_stars} stars...")

    # Step 1: Generate random star positions in a 2D unit square.
    # The uniform distribution and scale invariance mean the box size doesn't matter.
    stars = np.random.rand(num_stars, 2)

    # Step 2: Find the nearest neighbor for each star.
    # We use a k-d tree for efficient nearest neighbor searching.
    # We query for k=2 because the closest point (at distance 0) is the point itself.
    # The second closest is the actual nearest neighbor.
    kdtree = cKDTree(stars)
    # The `workers=-1` argument can speed up the query on multi-core systems.
    try:
        distances, indices = kdtree.query(stars, k=2, workers=-1)
    except TypeError: # Older scipy versions might not have `workers`
        distances, indices = kdtree.query(stars, k=2)

    # The 'successors' array maps each star to its nearest neighbor.
    # successors[i] is the index of the nearest neighbor to stars[i].
    successors = indices[:, 1]

    # Step 3: Build the underlying undirected graph.
    # The constellations are connected components, so we ignore edge direction.
    # An adjacency list is a good way to represent the graph for traversal.
    adj = [[] for _ in range(num_stars)]
    for i in range(num_stars):
        j = successors[i]
        adj[i].append(j)
        adj[j].append(i)

    # Step 4: Count the connected components (constellations).
    # We can use BFS or DFS. Here, we use BFS.
    visited = [False] * num_stars
    num_components = 0
    for i in range(num_stars):
        # If we haven't visited this star yet, it's part of a new component.
        if not visited[i]:
            num_components += 1
            # Start a BFS from this star to find all stars in its component.
            q = deque([i])
            visited[i] = True
            while q:
                u = q.popleft()
                for v in adj[u]:
                    if not visited[v]:
                        visited[v] = True
                        q.append(v)

    # Step 5: Calculate and print the final result.
    if num_components > 0:
        average_size = num_stars / num_components
        print(f"\nFound {num_components} constellations for {num_stars} stars.")
        print(f"The average number of stars per constellation is:")
        print(f"{num_stars} / {num_components} = {average_size}")
    else:
        print("No components found.")

    # For an infinite number of stars, this value converges to a known constant.
    # The known number of components is C = c_2 * N, where c_2 ≈ 0.3043.
    # The average size is N/C = 1/c_2 ≈ 3.286.
    print("\nThis is a simulation result. The true theoretical value for an infinite number of stars is ~3.286.")


solve_constellation_problem()