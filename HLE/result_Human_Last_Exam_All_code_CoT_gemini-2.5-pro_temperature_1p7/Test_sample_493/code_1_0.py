import numpy as np
import sys

# A simple Union-Find data structure to find connected components.
class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.num_sets = n

    def find(self, i):
        if self.parent[i] == i:
            return i
        self.parent[i] = self.find(self.parent[i])
        return self.parent[i]

    def union(self, i, j):
        root_i = self.find(i)
        root_j = self.find(j)
        if root_i != root_j:
            self.parent[root_i] = root_j
            self.num_sets -= 1

def distance_sq_torus(p1, p2):
    """Calculates squared Euclidean distance on a torus [0,1]x[0,1]."""
    dx = abs(p1[0] - p2[0])
    dy = abs(p1[1] - p2[1])
    # Account for the wrap-around nature of the torus
    dx = min(dx, 1 - dx)
    dy = min(dy, 1 - dy)
    return dx**2 + dy**2

def run_simulation(N):
    """
    Runs a single simulation for N stars and returns the average
    constellation size and the number of constellations found.
    """
    if N < 2:
        return N, 1 if N==1 else 0

    # 1. Generate N random stars on a 2D torus
    stars = np.random.rand(N, 2)

    # 2. Find the nearest neighbor for each star
    nearest_neighbors = np.empty(N, dtype=int)
    for i in range(N):
        min_dist_sq = float('inf')
        nn_index = -1
        for j in range(N):
            if i == j:
                continue
            dist_sq = distance_sq_torus(stars[i], stars[j])
            if dist_sq < min_dist_sq:
                min_dist_sq = dist_sq
                nn_index = j
        nearest_neighbors[i] = nn_index

    # 3. Identify connected components (constellations) using Union-Find
    uf = UnionFind(N)
    for i in range(N):
        uf.union(i, nearest_neighbors[i])

    # 4. Calculate average size for this simulation
    num_components = uf.num_sets
    if num_components == 0:
        return 0, 0
    
    avg_size = N / num_components
    return avg_size, num_components


if __name__ == '__main__':
    # Simulation parameters
    # Using a large N gives a more accurate result for the limit N -> infinity.
    # 2000 is a good balance of accuracy and speed for a demo.
    N = 2000 
    # More simulations lead to a more stable average.
    num_simulations = 10 

    all_results = []
    
    print(f"Running {num_simulations} simulations with {N} stars each.")
    print("This may take a minute...")
    
    for i in range(num_simulations):
        # We store the results of the last run to show in the final equation.
        last_avg_size, last_num_components = run_simulation(N)
        all_results.append(last_avg_size)
        # Progress indicator
        print(f"  Simulation {i+1}/{num_simulations} complete. Average size: {last_avg_size:.4f}")
        
    final_average = np.mean(all_results)
    
    print("\n-------------------------------------------")
    print("Example calculation from the final simulation:")
    # Using the numbers from the very last run as the example
    print(f"Number of Stars (N) = {N}")
    print(f"Number of Constellations = {last_num_components}")
    print(f"Final Equation: {N} / {last_num_components} = {last_avg_size:.4f}")
    print("-------------------------------------------")
    
    print(f"\nOverall average across all {num_simulations} simulations:")
    print(f"Average number of stars per constellation = {final_average:.4f}")
