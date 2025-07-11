import numpy as np
from scipy.spatial import cKDTree
import sys

# A higher recursion limit is set for the Union-Find algorithm to handle
# potentially deep call stacks with a large number of stars.
sys.setrecursionlimit(200000)

def solve_constellation_problem():
    """
    Simulates the constellation model to find the average number of stars
    per constellation.
    """
    # Number of stars for the simulation. A large number provides better statistics.
    N = 100000

    # A fixed random seed ensures that the simulation is reproducible.
    np.random.seed(42)

    # --- Union-Find Data Structure ---
    # This data structure is used to keep track of the connected components (constellations).
    parent = np.arange(N)
    size = np.ones(N, dtype=int)

    def find_set(i):
        """Finds the representative (root) of the set containing element i, with path compression."""
        if parent[i] == i:
            return i
        parent[i] = find_set(parent[i])
        return parent[i]

    def union_sets(i, j):
        """Merges the sets containing elements i and j, using union by size."""
        root_i = find_set(i)
        root_j = find_set(j)
        if root_i != root_j:
            # Attach the smaller tree to the root of the larger tree.
            if size[root_i] < size[root_j]:
                root_i, root_j = root_j, root_i
            parent[root_j] = root_i
            size[root_i] += size[root_j]

    # --- Simulation Steps ---

    # 1. Generate N stars uniformly in a 2D unit square.
    stars = np.random.rand(N, 2)

    # 2. Find the nearest neighbor for each star using a k-d tree.
    # We use a cKDTree with toroidal (periodic) boundary conditions to avoid edge effects.
    # The boxsize=[1, 1] argument tells the tree our space is a 1x1 square that wraps around.
    tree = cKDTree(stars, boxsize=[1, 1])

    # We query for k=2 nearest neighbors because the closest point to any star is the star itself.
    # The second closest point (at index 1) is its true nearest neighbor.
    distances, indices = tree.query(stars, k=2)
    nearest_neighbors = indices[:, 1]

    # 3. Build constellations by connecting each star to its nearest neighbor.
    # The union_sets function merges the components of a star and its neighbor.
    for i in range(N):
        union_sets(i, nearest_neighbors[i])

    # 4. Count the number of distinct constellations.
    # This is equal to the number of roots in the Union-Find structure, which are
    # the elements that are their own parent.
    num_constellations = 0
    for i in range(N):
        if parent[i] == i:
            num_constellations += 1

    # 5. Calculate and print the final result.
    average_size = N / num_constellations

    print("--- Simulation Results ---")
    print("The final equation for the average constellation size is: Total Stars / Number of Constellations")
    print(f"Total number of stars (N)       = {N}")
    print(f"Number of constellations found (C) = {num_constellations}")
    print(f"Average stars per constellation   = {N} / {num_constellations} = {average_size:.4f}")

if __name__ == '__main__':
    solve_constellation_problem()
