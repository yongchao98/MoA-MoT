import numpy as np
from scipy.spatial import KDTree

def solve_constellation_problem():
    """
    This function simulates the constellation model to find the average number
    of stars per constellation.
    """
    # Set a seed for reproducibility of the random point generation
    np.random.seed(42)

    # Number of stars (points) in the simulation.
    # A large number is used to get a result close to the theoretical average
    # for an infinite number of stars and to minimize boundary effects.
    N = 100000

    # 1. Generate N stars uniformly distributed in a 1x1 square.
    points = np.random.rand(N, 2)

    # 2. Use a KD-tree for efficient nearest neighbor search.
    # This is much faster than a brute-force O(N^2) comparison for large N.
    tree = KDTree(points)

    # For each point, find its nearest neighbor.
    # We query for k=2 because the closest point to any point is itself.
    # The second closest is its nearest neighbor in the dataset.
    distances, indices = tree.query(points, k=2)

    # The indices of the nearest neighbors are in the second column.
    nearest_neighbors = indices[:, 1]

    # 3. & 4. Find connected components using a Disjoint Set Union (DSU) data structure.
    class DSU:
        def __init__(self, n):
            self.parent = list(range(n))
            self._num_sets = n

        def find(self, i):
            if self.parent[i] == i:
                return i
            # Path compression for optimization
            self.parent[i] = self.find(self.parent[i])
            return self.parent[i]

        def union(self, i, j):
            root_i = self.find(i)
            root_j = self.find(j)
            if root_i != root_j:
                self.parent[root_i] = root_j
                self._num_sets -= 1

        def num_sets(self):
            return self._num_sets

    dsu = DSU(N)

    # For each star, form a connection (union) with its nearest neighbor.
    # This builds the graph structure defined by the model.
    for i in range(N):
        dsu.union(i, nearest_neighbors[i])

    # 4. The number of constellations is the number of disjoint sets.
    num_components = dsu.num_sets()

    # 5. Calculate the average number of stars per constellation.
    avg_size = N / num_components

    # Print the final equation with each number as requested.
    # Format: Total Stars / Number of Constellations = Average Size
    print(f"{N} / {num_components} = {avg_size}")

solve_constellation_problem()