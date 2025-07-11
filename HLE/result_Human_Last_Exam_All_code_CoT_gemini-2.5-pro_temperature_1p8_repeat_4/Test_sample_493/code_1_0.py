import numpy as np
from scipy.spatial import KDTree

def solve_constellation_problem():
    """
    Solves for the average number of stars per constellation via simulation.
    The code implements the plan described above.
    """
    # To ensure the simulation is reproducible, we set a fixed random seed.
    np.random.seed(42)

    # DSU (Disjoint Set Union) data structure for finding connected components.
    class DSU:
        """A simple and efficient Disjoint Set Union data structure."""
        def __init__(self, n):
            self.parent = list(range(n))
            self.num_sets = n

        def find(self, i):
            """Finds the representative of the set containing element i, with path compression."""
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

    # --- Simulation Parameters ---
    # Number of stars per simulation. Larger is better but slower. 20,000 is a good balance.
    N_STARS = 20000
    # Number of simulations to average over for a more stable result.
    NUM_SIMULATIONS = 10

    total_stars_accumulated = 0
    total_components_accumulated = 0

    print(f"Starting simulation...")
    print(f"Configuration: {NUM_SIMULATIONS} simulations of {N_STARS} stars each.")

    for i in range(NUM_SIMULATIONS):
        # 1. Generate N_STARS random 2D points in a unit square.
        points = np.random.rand(N_STARS, 2)

        # 2. Find the nearest neighbor for each point.
        # We use a KDTree with periodic boundary conditions (boxsize=[1,1])
        # to avoid edge effects. `query(points, k=2)` finds the 2 nearest
        # neighbors; the first is always the point itself.
        tree = KDTree(points, boxsize=[1, 1])
        # Use workers=-1 to leverage multiple CPU cores for speed.
        _, nearest_neighbor_indices = tree.query(points, k=2, workers=-1)
        
        # The actual nearest neighbor is the second one returned (index 1).
        neighbors = nearest_neighbor_indices[:, 1]

        # 3. Use DSU to count weakly connected components.
        # If star 'i' connects to its nearest neighbor 'j', we union their sets.
        dsu = DSU(N_STARS)
        for star_index in range(N_STARS):
            neighbor_index = neighbors[star_index]
            dsu.union(star_index, neighbor_index)

        # Accumulate results for the final average.
        total_stars_accumulated += N_STARS
        total_components_accumulated += dsu.num_sets
    
    print("Simulations complete.\n")

    # 4. Calculate the final average size based on all simulations.
    average_size = total_stars_accumulated / total_components_accumulated

    print("--- Final Calculation ---")
    print(f"Total stars simulated: {total_stars_accumulated}")
    print(f"Total constellations found: {total_components_accumulated}")
    print("\nEquation for the average constellation size:")
    # Here we output each number in the final equation as requested.
    print(f"{total_stars_accumulated} (stars) / {total_components_accumulated} (constellations) = {average_size:.4f}")
    
# Run the solver function
solve_constellation_problem()