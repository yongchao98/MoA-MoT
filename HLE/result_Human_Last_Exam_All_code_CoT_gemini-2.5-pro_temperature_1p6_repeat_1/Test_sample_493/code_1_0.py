import numpy as np
from scipy.spatial import KDTree

def solve_constellation_problem():
    """
    Solves the constellation problem via Monte Carlo simulation.
    """
    # Parameters for the simulation
    # N: Number of stars in each simulation
    # NUM_SIMULATIONS: Number of simulations to run and average
    N = 10000
    NUM_SIMULATIONS = 30

    print(f"Running {NUM_SIMULATIONS} simulations with {N} stars each to estimate the average constellation size...")

    all_avg_sizes = []
    last_run_details = {}

    for i in range(NUM_SIMULATIONS):
        # 1. Generate N random stars in a 2D space (e.g., a 1x1 square)
        points = np.random.rand(N, 2)

        # 2. Find the nearest neighbor for each star.
        # We use a KDTree for efficient searching (O(N log N)).
        # We query for k=2 because the closest point (k=1) to any point is itself.
        tree = KDTree(points)
        _, indices = tree.query(points, k=2)
        
        # adj[i] stores the index of the nearest neighbor of star i
        adj = indices[:, 1]

        # 3. Count the number of connected components (constellations).
        # A component is a set of stars connected by the nearest-neighbor relationship.
        # This is a directed graph where every node has an out-degree of 1.
        # The number of components is the number of cycles in this graph.
        visited = np.zeros(N, dtype=bool)
        num_components = 0
        for star_index in range(N):
            if not visited[star_index]:
                # Found the start of a new, unvisited component.
                num_components += 1
                # Traverse the component until we find a star we've already visited.
                current_star = star_index
                while not visited[current_star]:
                    visited[current_star] = True
                    current_star = adj[current_star]
        
        # 4. Calculate the average size for this simulation run.
        avg_size = N / num_components
        all_avg_sizes.append(avg_size)

        # Store the details of the last run to display later
        if i == NUM_SIMULATIONS - 1:
            last_run_details = {
                "N": N,
                "K": num_components,
                "avg_size": avg_size
            }

    # 5. Average the results from all simulations.
    final_average = np.mean(all_avg_sizes)

    print("\n--- Simulation Results ---")
    print(f"The estimated average number of stars per constellation is: {final_average:.4f}")
    
    # As requested, print the numbers for the final equation from the last run
    print("\nDetails from the final simulation run:")
    n_val = last_run_details['N']
    k_val = last_run_details['K']
    avg_val = last_run_details['avg_size']
    print(f"Total stars (N) = {n_val}")
    print(f"Number of constellations (K) = {k_val}")
    print(f"Final Equation: N / K = {n_val} / {k_val} = {avg_val:.4f}")

solve_constellation_problem()