import numpy as np
from scipy.spatial import cKDTree
from collections import deque

def solve_constellation_problem():
    """
    Simulates the constellation model to find the average number of stars per constellation.
    
    The theoretical value for this problem in 2D is known to be approximately 1.296.
    This simulation provides a computational method to approximate that value.
    """
    # Use a large number of stars for a more accurate result.
    N_STARS = 50000
    
    # We simulate in a unit square. The actual size doesn't matter due to uniform distribution.
    BOX_SIZE = 1.0

    # Step 1: Generate random star positions in a 2D box.
    points = np.random.rand(N_STARS, 2) * BOX_SIZE

    # Step 2: Find the nearest neighbor for each star using a k-d tree.
    # The `boxsize` argument handles toroidal (wrap-around) boundary conditions to avoid edge effects.
    # We query for k=2 because the closest point (k=1) to any point is the point itself.
    tree = cKDTree(points, boxsize=[BOX_SIZE, BOX_SIZE])
    distances, indices = tree.query(points, k=2)

    # The indices of the nearest neighbors are in the second column.
    nn_indices = indices[:, 1]

    # Step 3: Build the graph as an adjacency list.
    # The graph is undirected: if star i's nearest neighbor is j, they are connected.
    adj_list = [[] for _ in range(N_STARS)]
    for i in range(N_STARS):
        j = nn_indices[i]
        adj_list[i].append(j)
        adj_list[j].append(i)

    # Step 4: Find connected components (constellations) using Breadth-First Search (BFS).
    visited = [False] * N_STARS
    num_components = 0
    for i in range(N_STARS):
        if not visited[i]:
            # This is the start of a new, unvisited component.
            num_components += 1
            q = deque([i])
            visited[i] = True
            while q:
                u = q.popleft()
                # Find all neighbors of u and add them to the queue if not visited.
                for v in adj_list[u]:
                    if not visited[v]:
                        visited[v] = True
                        q.append(v)

    # Step 5: Calculate and print the average number of stars per constellation.
    avg_stars_per_constellation = N_STARS / num_components

    print("Running simulation for the nearest-neighbor constellation model...")
    print(f"Total number of stars = {N_STARS}")
    print(f"Number of constellations found = {num_components}")
    print("---")
    print(f"Average number of stars per constellation = {N_STARS} / {num_components} = {avg_stars_per_constellation}")

if __name__ == '__main__':
    solve_constellation_problem()
