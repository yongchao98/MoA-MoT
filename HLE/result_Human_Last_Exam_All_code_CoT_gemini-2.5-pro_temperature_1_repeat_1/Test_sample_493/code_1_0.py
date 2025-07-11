import numpy as np
from scipy.spatial import cKDTree

# This script requires the numpy and scipy libraries.
# You can install them using: pip install numpy scipy

class DSU:
    """
    A Disjoint Set Union (DSU) or Union-Find data structure.
    Used to efficiently track the connected components of a graph.
    """
    def __init__(self, n):
        # Initialize each node as its own parent, so we start with n sets.
        self.parent = list(range(n))
        self.num_sets = n

    def find(self, i):
        """Find the representative (root) of the set containing element i, with path compression."""
        if self.parent[i] == i:
            return i
        # Path compression: set the parent of i directly to the root.
        self.parent[i] = self.find(self.parent[i])
        return self.parent[i]

    def union(self, i, j):
        """Merge the sets containing elements i and j."""
        root_i = self.find(i)
        root_j = self.find(j)
        if root_i != root_j:
            # Union by making one root the parent of the other.
            self.parent[root_i] = root_j
            self.num_sets -= 1
            return True
        return False

def calculate_average_constellation_size(num_stars, box_size):
    """
    Simulates the constellation model and calculates the average size.
    
    Args:
        num_stars (int): The total number of stars to simulate.
        box_size (float): The side length of the square area for the simulation.

    Returns:
        A tuple containing (num_stars, num_components, average_size).
    """
    # 1. Generate star positions
    # We generate 'num_stars' points uniformly in a square of size 'box_size' x 'box_size'.
    # This simulates the "uniformly distributed in the sky" condition.
    # Using a large number of stars helps minimize edge effects and get closer to the
    # theoretical asymptotic value.
    points = np.random.rand(num_stars, 2) * box_size

    # 2. Find the nearest neighbor for each star
    # We use scipy's cKDTree for a very fast nearest-neighbor search.
    # We query for the 2 nearest neighbors because the first nearest neighbor (k=1)
    # to any point is always the point itself. The second (k=2) is the one we want.
    tree = cKDTree(points)
    # `indices` will be an array where `indices[i, 1]` is the index of the
    # nearest neighbor to the point `points[i]`.
    distances, indices = tree.query(points, k=2, workers=-1)
    nn_indices = indices[:, 1]

    # 3. Count the number of constellations (connected components)
    # We model this as a graph problem. Each star is a vertex.
    # An undirected edge exists between two stars if one is the nearest neighbor of the other.
    # A constellation is a connected component in this graph.
    # The DSU data structure is perfect for finding the number of components.
    dsu = DSU(num_stars)
    for i in range(num_stars):
        neighbor_index = nn_indices[i]
        # We merge the set containing star 'i' and the set containing its nearest neighbor.
        dsu.union(i, neighbor_index)

    num_components = dsu.num_sets
    
    # 4. Calculate the average number of stars per constellation
    if num_components == 0:
        return num_stars, 0, 0
    average_size = num_stars / num_components
    
    return num_stars, num_components, average_size

# --- Main execution ---
if __name__ == "__main__":
    # A large number of stars is used to ensure the result is close to the theoretical
    # value for an infinite number of stars. The box size is chosen to keep the
    # star density constant at 1.0, though the final result is independent of density.
    num_stars_to_simulate = 50000
    simulation_area_size = np.sqrt(num_stars_to_simulate)

    N, K, avg_size = calculate_average_constellation_size(
        num_stars_to_simulate, 
        simulation_area_size
    )

    # Output the results of the simulation.
    # The prompt asks to "output each number in the final equation".
    print(f"Simulation Parameters:")
    print(f"Total number of stars (N): {N}")
    print(f"Number of constellations found (K): {K}")
    print("\nFinal Equation:")
    print(f"Average stars per constellation = N / K")
    print(f"{N} / {K} = {avg_size:.4f}")