import numpy as np
from scipy.spatial import KDTree

class DSU:
    """
    A simple Disjoint Set Union (DSU) data structure with path compression.
    Used to efficiently count connected components in a graph.
    """
    def __init__(self, n):
        self.parent = list(range(n))
        self.num_sets = n

    def find(self, i):
        """Find the representative (root) of the set containing element i."""
        if self.parent[i] == i:
            return i
        # Path compression
        self.parent[i] = self.find(self.parent[i])
        return self.parent[i]

    def union(self, i, j):
        """Merge the sets containing elements i and j."""
        root_i = self.find(i)
        root_j = self.find(j)
        if root_i != root_j:
            # A simple merge, could be improved with union-by-rank/size
            self.parent[root_i] = root_j
            self.num_sets -= 1

def run_constellation_simulation(num_stars, seed=42):
    """
    Simulates the constellation model and calculates the average size.
    
    Args:
        num_stars (int): The number of stars to simulate.
        seed (int): A seed for the random number generator for reproducibility.
        
    Returns:
        float: The average number of stars per constellation.
    """
    # Set seed for reproducible results
    np.random.seed(seed)
    
    # 1. Generate points on a 2D torus [0,1) x [0,1)
    # The boxsize argument to KDTree handles the toroidal (periodic) boundary conditions.
    points = np.random.rand(num_stars, 2)
    
    # 2. Find the nearest neighbor for each point
    # We use a KDTree for efficient nearest neighbor search in O(N log N).
    # The 'boxsize' argument tells the KDTree to use toroidal distances.
    kdtree = KDTree(points, boxsize=[1.0, 1.0])
    
    # Query for the 2 nearest neighbors because the closest point to any point
    # is the point itself. The second closest is the one we want.
    _, indices = kdtree.query(points, k=2)
    
    # The nearest neighbor index for point i is in indices[i, 1]
    nearest_neighbors = indices[:, 1]
    
    # 3. Use DSU to find connected components (constellations)
    dsu = DSU(num_stars)
    for i in range(num_stars):
        neighbor_index = nearest_neighbors[i]
        dsu.union(i, neighbor_index)
        
    # 4. Calculate the average size of a constellation
    num_constellations = dsu.num_sets
    average_size = num_stars / num_constellations
    
    print("--- Simulation Results ---")
    print(f"Total number of stars (N): {num_stars}")
    print(f"Number of constellations found (k): {num_constellations}")
    print("\nFinal Calculation:")
    print(f"Average constellation size = N / k")
    print(f"{num_stars} / {num_constellations} = {average_size:.4f}")
    
    return average_size

if __name__ == '__main__':
    # Number of stars for the simulation.
    # A larger number provides a more accurate estimate.
    N = 50000
    
    run_constellation_simulation(N)
