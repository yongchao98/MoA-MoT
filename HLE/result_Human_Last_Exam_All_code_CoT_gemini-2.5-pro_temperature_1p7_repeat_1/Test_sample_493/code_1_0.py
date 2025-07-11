import numpy as np
from scipy.spatial import cKDTree

def solve_constellation_problem():
    """
    Simulates the constellation model to find the average number of stars
    per constellation.
    """

    class UnionFind:
        """
        An efficient implementation of the Union-Find data structure
        used to count connected components in a graph.
        """
        def __init__(self, n):
            self.parent = list(range(n))
            self.num_sets = n

        def find(self, i):
            """Finds the representative of the set containing element i with path compression."""
            if self.parent[i] == i:
                return i
            self.parent[i] = self.find(self.parent[i])
            return self.parent[i]

        def union(self, i, j):
            """Merges the sets containing elements i and j."""
            root_i = self.find(i)
            root_j = self.find(j)
            if root_i != root_j:
                self.parent[root_i] = root_j
                self.num_sets -= 1
                return True
            return False

    # Simulation Parameters
    # A larger number of stars leads to a more accurate result.
    N_STARS = 200000
    # The size of the simulation area (a square). We use a unit box.
    BOX_SIZE = 1.0
    # Use a fixed seed for reproducibility.
    SEED = 42
    
    print(f"Running simulation with {N_STARS} stars...")
    np.random.seed(SEED)

    # 1. Generate N random stars in a 2D square.
    # We use a uniform distribution on [0, BOX_SIZE) x [0, BOX_SIZE).
    stars = np.random.uniform(0, BOX_SIZE, size=(N_STARS, 2))

    # 2. Find the nearest neighbor for each star.
    # We build a KD-Tree with periodic boundary conditions (boxsize)
    # to simulate an infinite, uniform distribution of stars without edge effects.
    # We query for k=2 because the closest point (k=1) to any star is itself.
    print("- Finding nearest neighbor for each star...")
    tree = cKDTree(stars, boxsize=[BOX_SIZE, BOX_SIZE])
    _, indices = tree.query(stars, k=2)
    nearest_neighbors = indices[:, 1]

    # 3. Count constellations using the Union-Find data structure.
    # A constellation is a connected component in the undirected graph
    # where an edge exists between each star and its nearest neighbor.
    print("- Identifying constellations (connected components)...")
    uf = UnionFind(N_STARS)
    for i in range(N_STARS):
        uf.union(i, nearest_neighbors[i])

    num_constellations = uf.num_sets
    
    # 4. Calculate the average size.
    avg_size = N_STARS / num_constellations

    print("\n--- Simulation Results ---")
    print(f"Total number of stars simulated: {N_STARS}")
    print(f"Number of constellations found: {num_constellations}")
    print("\nAverage number of stars per constellation is given by:")
    print(f"Total Stars / Number of Constellations")
    print(f"Final equation: {N_STARS} / {num_constellations} = {avg_size}")


if __name__ == '__main__':
    solve_constellation_problem()