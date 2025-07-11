import numpy as np
from scipy.spatial import cKDTree
from collections import deque
import sys

def run_simulation(num_stars, num_trials):
    """
    Runs a Monte Carlo simulation to find the average size of a constellation
    in a random nearest-neighbor graph.
    """
    all_trial_results = []

    print(f"Running {num_trials} trials with {num_stars} stars each...")
    # A simple text-based progress bar
    sys.stdout.write("Progress: [%s]" % (" " * num_trials))
    sys.stdout.flush()
    sys.stdout.write("\b" * (num_trials + 1)) # return to start of line, after '['

    # Store results from the last trial to display in the equation
    last_trial_components = 0

    for trial in range(num_trials):
        # 1. Generate N uniformly distributed points in a 2D square
        points = np.random.rand(num_stars, 2)

        # 2. Find the nearest neighbor for each point using a KD-Tree for efficiency.
        # This creates a directed graph where each node has an out-degree of 1.
        tree = cKDTree(points)
        # Query for k=2 because the closest point (k=1) to any point is itself.
        # We want the second-closest point.
        distances, indices = tree.query(points, k=2)
        
        # The nearest neighbor of point `i` is `indices[i, 1]`
        nearest_neighbors = indices[:, 1]

        # 3. Build an undirected representation of the graph to find connected components.
        # Two stars are in the same constellation if there's any path between them.
        undirected_adj = {i: [] for i in range(num_stars)}
        for i in range(num_stars):
            j = nearest_neighbors[i]
            # Add an edge in both directions to allow traversal across the component
            undirected_adj[i].append(j)
            undirected_adj[j].append(i)

        # 4. Count the number of weakly connected components using BFS.
        # Each component corresponds to one constellation.
        visited = set()
        num_components = 0
        for i in range(num_stars):
            if i not in visited:
                num_components += 1
                q = deque([i])
                visited.add(i)
                while q:
                    u = q.popleft()
                    for v in undirected_adj[u]:
                        if v not in visited:
                            visited.add(v)
                            q.append(v)
        
        last_trial_components = num_components
        
        # 5. Calculate the average constellation size for this trial and store it.
        avg_size_this_trial = num_stars / num_components
        all_trial_results.append(avg_size_this_trial)

        # Update progress bar
        sys.stdout.write("=")
        sys.stdout.flush()

    sys.stdout.write("]\n") # End of progress bar

    # Calculate the final average over all trials for a stable estimate
    final_average = np.mean(all_trial_results)

    return final_average, num_stars, last_trial_components

if __name__ == '__main__':
    # --- Parameters for the simulation ---
    # A larger number of stars reduces edge effects and better approximates an infinite plane.
    NUM_STARS = 5000
    # A larger number of trials provides a more stable and reliable average.
    NUM_TRIALS = 50

    # Run the simulation
    final_avg, stars_for_eq, comps_for_eq = run_simulation(NUM_STARS, NUM_TRIALS)

    # --- Print the results ---
    print("\n--- Simulation Results ---")
    print(f"The model connects each star to its single nearest neighbor.")
    print(f"A 'constellation' is a resulting group of connected stars.")
    
    print("\nBased on the final trial of the simulation:")
    print("The average number of stars per constellation is calculated as:")
    print(f"Average Size = Total Stars / Number of Constellations")
    print(f"             = {stars_for_eq} / {comps_for_eq}")
    print(f"             = {stars_for_eq / comps_for_eq:.4f}")

    print("\nAfter averaging the results over all trials, the estimated value is:")
    print(f"Final Average Stars per Constellation: {final_avg:.4f}")
    
    print("\nNote: The simulation provides an estimate. The problem has a known theoretical")
    print("value in stochastic geometry, which this simulation approximates.")
