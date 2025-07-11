import numpy as np
from scipy.spatial import KDTree
import sys

def solve():
    """
    Performs a Monte Carlo simulation to find the average number of stars
    per constellation in a nearest-neighbor graph model.
    """
    # It's good practice to increase the recursion limit for Union-Find,
    # as it can involve deep recursive calls for large datasets.
    sys.setrecursionlimit(30000)

    def count_components_for_trial(num_points):
        """
        Generates N points, finds nearest neighbors, and counts graph components.
        
        Args:
            num_points (int): The number of stars to simulate in this trial.

        Returns:
            int: The number of connected components (constellations).
        """
        # 1. Generate N random points (stars) in a [0, 1] x [0, 1] square.
        points = np.random.rand(num_points, 2)

        # 2. Use a KDTree to efficiently find the nearest neighbor for each point.
        # We query for k=2 because the closest point to any point is itself.
        tree = KDTree(points)
        _, indices = tree.query(points, k=2)
        
        # The nearest neighbor of point `i` is the second column of the result.
        nearest_neighbors = indices[:, 1]

        # 3. Count the connected components using the Union-Find algorithm.
        parent = list(range(num_points))
        num_components = num_points

        def find_set(v):
            # Find the root of the set containing element v with path compression.
            if v == parent[v]:
                return v
            parent[v] = find_set(parent[v])
            return parent[v]

        # Iterate through all stars, connecting each to its nearest neighbor.
        for i in range(num_points):
            # Find the sets of the star and its neighbor.
            set_i = find_set(i)
            set_j = find_set(nearest_neighbors[i])
            
            # If they are not already in the same set, merge them.
            if set_i != set_j:
                parent[set_j] = set_i
                num_components -= 1
        
        return num_components

    # --- Simulation Parameters ---
    # Number of stars in each individual simulation. A larger number
    # reduces boundary effects and better approximates the infinite plane.
    num_points_per_trial = 10000
    
    # Number of trials to run. Averaging over many trials gives a more stable result.
    num_trials = 100

    print("Running Monte Carlo simulation to estimate the average constellation size...")
    print(f"Configuration: {num_trials} trials, each with {num_points_per_trial} stars.")
    print("-" * 20)

    total_stars_simulated = 0
    total_components_found = 0

    # --- Main Simulation Loop ---
    for i in range(num_trials):
        total_stars_simulated += num_points_per_trial
        
        components_in_trial = count_components_for_trial(num_points_per_trial)
        total_components_found += components_in_trial
        
        # Print progress to keep the user informed.
        if (i + 1) % 10 == 0:
            print(f"Completed trial {i+1}/{num_trials}...")

    # --- Final Calculation ---
    # The average constellation size is the total number of stars observed
    # across all trials divided by the total number of constellations found.
    average_size = total_stars_simulated / total_components_found
    
    print("\n--- Simulation Complete ---")
    print("The final estimate is based on the equation:")
    print("Average Size = Total Stars Simulated / Total Constellations Found")
    print("\nFinal Equation Details:")
    print(f"Total Stars Simulated   = {total_stars_simulated}")
    print(f"Total Constellations Found = {total_components_found}")
    print(f"Average Stars per Constellation = {total_stars_simulated} / {total_components_found} = {average_size:.4f}")

# Execute the main function
solve()