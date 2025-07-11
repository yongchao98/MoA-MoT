import numpy as np
from scipy.spatial import KDTree

def solve_constellation_problem():
    """
    This function simulates the star constellation problem to find the average
    number of stars per constellation.
    """
    # Step 1: Generate a large number of stars (points)
    # A larger N gives a more accurate approximation of the true average.
    # We use a fixed random seed for reproducibility of the simulation result.
    N = 50000
    np.random.seed(42)
    points = np.random.rand(N, 2)

    # Step 2: For each star, find its nearest neighbor
    # We use a KDTree data structure for an efficient search, which is much
    # faster than comparing every point with every other point.
    # tree.query(points, k=2) finds the two nearest neighbors for each point.
    # The first (at index 0) is the point itself, so its nearest neighbor is
    # the second one (at index 1).
    tree = KDTree(points)
    _, nn_indices = tree.query(points, k=2)
    nearest_neighbor_indices = nn_indices[:, 1]

    # Step 3: Build the graph representation (adjacency list)
    # The constellations are connected components in an undirected graph where
    # each star is joined to its nearest neighbor.
    adj = [[] for _ in range(N)]
    for i in range(N):
        j = nearest_neighbor_indices[i]
        # Add an undirected edge between star i and its nearest neighbor j
        adj[i].append(j)
        adj[j].append(i)

    # Step 4: Count the connected components (constellations)
    # We use a traversal algorithm (an iterative DFS) to find the components.
    visited = [False] * N
    num_components = 0
    for i in range(N):
        if not visited[i]:
            num_components += 1
            stack = [i]
            visited[i] = True
            while stack:
                u = stack.pop()
                for v in adj[u]:
                    if not visited[v]:
                        visited[v] = True
                        stack.append(v)

    # Step 5: Calculate and print the average size of a constellation
    average_size = N / num_components

    print(f"Simulation run with N = {N} stars.")
    print(f"Number of constellations found = {num_components}.")
    print("\nFinal Equation:")
    print(f"{N} / {num_components} = {average_size}")


if __name__ == '__main__':
    solve_constellation_problem()

<<<3.248473232848233>>>