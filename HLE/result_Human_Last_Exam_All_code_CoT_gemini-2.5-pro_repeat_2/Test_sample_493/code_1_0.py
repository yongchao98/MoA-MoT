import numpy as np
from scipy.spatial import cKDTree
from collections import deque

def simulate_constellations():
    """
    This script simulates the formation of "constellations" based on a nearest-neighbor
    model and calculates the average number of stars per constellation.
    """
    # --- Simulation Parameters ---
    # Number of stars to generate in each simulation.
    # A larger number gives a more accurate result but is slower.
    N_STARS = 20000

    # Number of times to run the simulation to average the results.
    # This improves the stability and accuracy of the final answer.
    N_TRIALS = 10

    print(f"Starting simulation with {N_STARS} stars and {N_TRIALS} trials.")
    print("-" * 30)

    total_avg_size = 0.0
    last_trial_details = {}

    for trial in range(N_TRIALS):
        # 1. Generate N_STARS random points in a 2D unit square.
        points = np.random.rand(N_STARS, 2)

        # 2. Use a cKDTree for efficient nearest-neighbor lookup.
        # We query for k=2 because the nearest point to any point is itself.
        tree = cKDTree(points)
        _, indices = tree.query(points, k=2)
        
        # The second column contains the index of the nearest neighbor for each point.
        nearest_indices = indices[:, 1]

        # 3. Build the adjacency list for the UNDIRECTED graph.
        # An edge exists between u and v if v=nn(u) or u=nn(v).
        adj = [set() for _ in range(N_STARS)]
        for i in range(N_STARS):
            j = nearest_indices[i]
            # Add the undirected edge {i, j}
            adj[i].add(j)
            adj[j].add(i)

        # 4. Count the connected components (constellations) using BFS.
        visited = [False] * N_STARS
        num_components = 0
        for i in range(N_STARS):
            if not visited[i]:
                num_components += 1
                q = deque([i])
                visited[i] = True
                while q:
                    u = q.popleft()
                    for v in adj[u]:
                        if not visited[v]:
                            visited[v] = True
                            q.append(v)
        
        # 5. Calculate the average size for this trial.
        if num_components > 0:
            avg_size_this_trial = N_STARS / num_components
            total_avg_size += avg_size_this_trial
            print(f"Trial {trial + 1}/{N_TRIALS}: Found {num_components} constellations. Average size = {avg_size_this_trial:.4f}")
            if trial == N_TRIALS - 1:
                last_trial_details = {
                    "stars": N_STARS,
                    "components": num_components,
                    "avg_size": avg_size_this_trial
                }

    # 6. Calculate and print the final average over all trials.
    final_avg_size = total_avg_size / N_TRIALS

    print("-" * 30)
    print("Final Calculation Details (from last trial):")
    n = last_trial_details['stars']
    c = last_trial_details['components']
    a = last_trial_details['avg_size']
    print(f"Average Stars per Constellation: {n} / {c} = {a:.4f}")
    print("-" * 30)
    print(f"Stable Average Size (over {N_TRIALS} trials): {final_avg_size:.4f}")


if __name__ == "__main__":
    simulate_constellations()