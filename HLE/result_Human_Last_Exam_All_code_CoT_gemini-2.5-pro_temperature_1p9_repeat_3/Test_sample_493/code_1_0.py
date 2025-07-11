import numpy as np
from sklearn.neighbors import NearestNeighbors

def solve_constellation_problem():
    """
    Simulates the star constellation model to find the average number of stars
    per constellation.
    """
    # --- Parameters ---
    # Number of stars in each simulation. A large number reduces boundary effects.
    N_STARS = 10000
    # Number of simulations to run. Averaging results gives a more stable estimate.
    NUM_SIMULATIONS = 30

    # --- Simulation Loop ---
    all_run_averages = []
    print(f"Starting simulation with N={N_STARS} stars, running {NUM_SIMULATIONS} times.")
    print("-" * 40)

    for i in range(NUM_SIMULATIONS):
        # 1. Generate N random stars (points) in a 2D plane (unit square)
        points = np.random.rand(N_STARS, 2)

        # 2. For each star, find its nearest neighbor
        # n_neighbors=2 is used because the closest point to any point is itself.
        # The nearest *neighbor* is the second point returned.
        nbrs = NearestNeighbors(n_neighbors=2, algorithm='auto').fit(points)
        _, indices = nbrs.kneighbors(points)
        
        # next_node[j] stores the index of the nearest neighbor of star j
        next_node = indices[:, 1]

        # 3. Find the number of constellations by counting cycles
        # We use a component_id array to keep track of which component each star belongs to.
        # A value of 0 means it's unassigned.
        component_id = np.zeros(N_STARS, dtype=int)
        num_components = 0

        for star_idx in range(N_STARS):
            if component_id[star_idx] == 0:
                # This star has not been assigned to a component yet.
                # We trace its path of nearest neighbors.
                
                path_dict = {}  # Using a dict for fast O(1) lookup of nodes in path
                path_list = []  # Keep track of path order
                current_star = star_idx
                
                # Follow the path until we find a star that is already in a component
                # or a star that is already on our current path (a new cycle).
                while current_star not in path_dict and component_id[current_star] == 0:
                    path_dict[current_star] = len(path_list)
                    path_list.append(current_star)
                    current_star = next_node[current_star]
                
                # The path trace has ended. Determine why.
                if current_star in path_dict:
                    # We found a new cycle because we looped back onto our current path.
                    # This defines a new component.
                    num_components += 1
                    id_to_assign = num_components
                    
                    # All stars on this path belong to the new component.
                    for node in path_list:
                        component_id[node] = id_to_assign
                else:
                    # We merged into a pre-existing component.
                    id_to_assign = component_id[current_star]

                    # All stars on this path lead to this existing component.
                    for node in path_list:
                        component_id[node] = id_to_assign
        
        # 4. Calculate the average size for this simulation run
        if num_components > 0:
            avg_size_this_run = N_STARS / num_components
            all_run_averages.append(avg_size_this_run)

            # As requested, output the equation for the last run
            if i == NUM_SIMULATIONS - 1:
                 print("\nSample calculation from the last run:")
                 print(f"{N_STARS} stars / {num_components} constellations = {avg_size_this_run:.4f} average size")
                 print("-" * 40)

    # 5. Calculate and print the final average over all simulations
    final_average = np.mean(all_run_averages)
    print("\nFinal Result:")
    print(f"The average number of stars per constellation is approximately: {final_average:.4f}")

if __name__ == '__main__':
    solve_constellation_problem()
