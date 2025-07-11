import numpy as np
from scipy.spatial import cKDTree

def estimate_average_constellation_size():
    """
    This function estimates the average number of stars per constellation
    by simulating the described model.
    """
    # --- Simulation Parameters ---
    # Number of stars in each simulation. A larger number gives a more accurate result.
    N_STARS = 5000
    # Number of simulations to run. More simulations lead to a more stable average.
    N_SIMULATIONS = 20
    
    # These will accumulate totals across all simulations
    total_stars_simulated = 0
    total_constellations_found = 0

    print(f"Running {N_SIMULATIONS} simulations with {N_STARS} stars each...")

    for i in range(N_SIMULATIONS):
        # 1. Generate N_STARS random points (stars) in a 2D unit square.
        points = np.random.rand(N_STARS, 2)
        
        # 2. Find the nearest neighbor for each point.
        # We use a cKDTree for efficient nearest neighbor searching.
        # k=2 because the first nearest neighbor of a point is the point itself.
        try:
            tree = cKDTree(points)
            # distances, indices of nearest neighbors
            _, nn_indices = tree.query(points, k=2) 
            # The nearest neighbor is the second column (index 1)
            nearest_neighbor_map = nn_indices[:, 1]
        except ImportError:
            print("SciPy is not installed. Please install it with 'pip install scipy'")
            return

        # 3. Build the undirected graph from the nearest neighbor connections.
        # We use an adjacency list represented by a dictionary.
        adj = {star: set() for star in range(N_STARS)}
        for star_idx, nn_idx in enumerate(nearest_neighbor_map):
            # Add an undirected edge between a star and its nearest neighbor
            adj[star_idx].add(nn_idx)
            adj[nn_idx].add(star_idx)

        # 4. Find the number of connected components (constellations) using BFS.
        visited = [False] * N_STARS
        num_constellations = 0
        
        for star_idx in range(N_STARS):
            if not visited[star_idx]:
                # This is the start of a new, unvisited constellation.
                num_constellations += 1
                q = [star_idx]  # A queue for BFS
                visited[star_idx] = True
                
                # Explore the entire component
                while q:
                    current_star = q.pop(0)
                    for neighbor in adj[current_star]:
                        if not visited[neighbor]:
                            visited[neighbor] = True
                            q.append(neighbor)
        
        # 5. Add this simulation's results to the totals.
        total_stars_simulated += N_STARS
        total_constellations_found += num_constellations
        
        print(f"  Simulation {i+1}/{N_SIMULATIONS} complete.")


    # 6. Calculate the final average and print the result.
    if total_constellations_found > 0:
        average_size = total_stars_simulated / total_constellations_found
        
        print("\n--- Simulation Results ---")
        print(f"Total stars simulated: {total_stars_simulated}")
        print(f"Total constellations found: {total_constellations_found}")
        print("\nAverage number of stars per constellation (Total Stars / Total Constellations):")
        
        # The final equation and its result
        print(f"{total_stars_simulated} / {total_constellations_found} = {average_size}")
        
        # Required final answer format
        print(f"\n<<<{average_size:.4f}>>>")
    else:
        print("Could not calculate the average size.")

if __name__ == '__main__':
    # NOTE: This script requires the 'numpy' and 'scipy' libraries.
    # You can install them using pip:
    # pip install numpy scipy
    estimate_average_constellation_size()