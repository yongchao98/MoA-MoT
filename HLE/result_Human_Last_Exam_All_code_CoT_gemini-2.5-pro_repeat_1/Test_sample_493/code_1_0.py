import numpy as np
from scipy.spatial import KDTree

# A simple Union-Find implementation for counting connected components
class UnionFind:
    """
    A data structure for efficiently tracking connected components in a graph.
    """
    def __init__(self, n):
        self.parent = list(range(n))
        self.num_sets = n

    def find(self, i):
        # Find the root of the set containing element i with path compression
        if self.parent[i] == i:
            return i
        self.parent[i] = self.find(self.parent[i])
        return self.parent[i]

    def union(self, i, j):
        # Merge the sets containing elements i and j
        root_i = self.find(i)
        root_j = self.find(j)
        if root_i != root_j:
            self.parent[root_j] = root_i
            self.num_sets -= 1
            return True
        return False

def run_simulation(n_points):
    """
    Runs a single simulation for a given number of points.
    Returns the average constellation size for this run.
    """
    # 1. Generate N random points in a 2D square
    points = np.random.rand(n_points, 2)

    # 2. Find the nearest neighbor for each point using a KD-Tree for efficiency.
    # We query for k=2 because the point itself is its own nearest neighbor at distance 0.
    tree = KDTree(points)
    _, nn_indices = tree.query(points, k=2)
    
    # The actual nearest neighbor's index is in the second column.
    nearest_neighbors = nn_indices[:, 1]

    # 3. Build the undirected graph and count connected components using Union-Find.
    # For each star and its nearest neighbor, we merge their sets.
    uf = UnionFind(n_points)
    for i in range(n_points):
        j = nearest_neighbors[i]
        uf.union(i, j)

    # 4. Calculate the average size of a component (constellation).
    num_components = uf.num_sets
    if num_components == 0:
        return 0
    return n_points / num_components

# --- Main execution ---

# Simulation Parameters
N_POINTS = 5000
N_SIMULATIONS = 20  # Lower for speed, higher for accuracy

print(f"Running simulation with {N_POINTS} stars over {N_SIMULATIONS} trials...")
all_results = [run_simulation(N_POINTS) for _ in range(N_SIMULATIONS)]
simulation_avg_size = np.mean(all_results)
print(f"Average stars per constellation (from simulation): {simulation_avg_size:.4f}")
print("-" * 20)

# Theoretical Calculation
print("Theoretical calculation based on asymptotic results:")
print("The average constellation size 'S' for a large number of stars is given by:")
print("S = 2 / p")
print("where 'p' is the probability that a star is a mutual nearest neighbor to its nearest neighbor.")
print("The value of 'p' for a 2D uniform distribution is p = 6*pi / (8*pi + sqrt(3)).")
print("\nLet's calculate the value step-by-step:")

# Breaking down the formula
six_pi = 6 * np.pi
eight_pi = 8 * np.pi
sqrt_3 = np.sqrt(3)

# Final probability p
p = six_pi / (eight_pi + sqrt_3)

print(f"p = (6 * {np.pi:.4f}) / (8 * {np.pi:.4f} + {sqrt_3:.4f})")
print(f"p = {six_pi:.4f} / ({eight_pi:.4f} + {sqrt_3:.4f})")
print(f"p = {six_pi:.4f} / {eight_pi + sqrt_3:.4f}")
print(f"p ≈ {p:.4f}")

# The average size is 2/p, which can be simplified algebraically to:
# S = 2 / p = 2 * (8*pi + sqrt(3)) / (6*pi) = (8*pi + sqrt(3)) / (3*pi) = 8/3 + sqrt(3)/(3*pi)
print("\nThe final formula for the average size 'S' is:")
print("S = 8/3 + sqrt(3) / (3 * pi)")

val1 = 8/3
val2 = sqrt_3
val3 = 3 * np.pi
theoretical_avg_size = val1 + val2 / val3

print(f"S = {val1:.4f} + {val2:.4f} / ({val3/np.pi:.0f} * {np.pi:.4f})")
print(f"S = {val1:.4f} + {val2:.4f} / {val3:.4f}")
print(f"S = {val1:.4f} + {val2/val3:.4f}")
print(f"S ≈ {theoretical_avg_size:.4f}")

<<<2.8505>>>