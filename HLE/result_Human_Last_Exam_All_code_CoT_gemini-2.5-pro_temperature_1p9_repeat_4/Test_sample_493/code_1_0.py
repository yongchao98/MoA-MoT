import numpy as np
from scipy.spatial import cKDTree
import collections

def simulate_and_calculate_average(num_stars, num_simulations):
    """
    Runs a Monte Carlo simulation to find the average size of a constellation.

    In this model, stars are points in a 2D plane, and each star is connected
    to its nearest neighbor, forming a graph. A constellation is a weakly
    connected component of this graph.

    Args:
        num_stars (int): The number of stars to use in each simulation.
                         A larger number gives a better approximation.
        num_simulations (int): The number of times to run the simulation
                               to average the results for a stable estimate.

    Returns:
        float: The estimated average number of stars per constellation.
    """
    total_avg_size = 0.0
    last_avg = 0.0
    num_components = 0

    print(f"Running {num_simulations} simulations with {num_stars} stars each...")

    for i in range(num_simulations):
        # 1. Generate star positions in a square with side length chosen
        # to keep the density constant. This simulates a uniform distribution.
        # We use periodic boundary conditions to avoid edge effects.
        box_size = np.sqrt(num_stars)
        points = np.random.rand(num_stars, 2) * box_size

        # 2. Find the nearest neighbor for each star using a cKDTree.
        # cKDTree is highly optimized for this and handles periodic boundaries
        # via the 'boxsize' parameter. We query for k=2 because the
        # closest point to any given point is itself.
        kdtree = cKDTree(points, boxsize=[box_size, box_size])
        _, indices = kdtree.query(points, k=2, workers=-1)
        # The nearest neighbor indices are in the second column.
        nearest_neighbors = indices[:, 1]

        # 3. Build an undirected graph to find the weakly connected components.
        # If A's nearest neighbor is B, A and B are in the same component.
        adj = collections.defaultdict(list)
        for star_idx, neighbor_idx in enumerate(nearest_neighbors):
            adj[star_idx].append(neighbor_idx)
            adj[neighbor_idx].append(star_idx)

        # 4. Use BFS to find all connected components ("constellations").
        visited = [False] * num_stars
        num_components = 0
        for star_idx in range(num_stars):
            if not visited[star_idx]:
                num_components += 1
                q = collections.deque([star_idx])
                visited[star_idx] = True
                while q:
                    u = q.popleft()
                    for v in adj[u]:
                        if not visited[v]:
                            visited[v] = True
                            q.append(v)
        
        # 5. Calculate the average size for this single simulation.
        if num_components > 0:
            last_avg = num_stars / num_components
            total_avg_size += last_avg
        
        # Optional: Print progress
        # print(f"  Simulation {i+1}/{num_simulations} -> Avg size: {last_avg:.4f}")

    # Calculate the final average over all simulations
    final_average = total_avg_size / num_simulations

    # Print the final breakdown as requested
    print("\n--- Final Results ---")
    print("This result is based on the final simulation run:")
    equation = f"Total Stars ({num_stars}) / Number of Constellations ({num_components}) = Average Size ({last_avg:.4f})"
    print(equation)
    print(f"\nFinal estimated average stars per constellation after {num_simulations} simulations: {final_average:.4f}")

    return final_average

if __name__ == '__main__':
    # Simulation parameters: More stars and more simulations yield a more accurate result.
    # These values provide a good balance between accuracy and runtime.
    NUM_STARS = 5000
    NUM_SIMULATIONS = 50
    
    # Run the simulation and get the final answer.
    average_constellation_size = simulate_and_calculate_average(NUM_STARS, NUM_SIMULATIONS)
    
    # The known theoretical value is ~3.217. The simulation should be close to this.
