import numpy as np

def solve_constellation_problem():
    """
    Simulates the star constellation problem to find the average number of stars
    per constellation.
    """
    # Simulation parameters
    # N: Number of stars in each simulation. A larger number gives a better
    #    approximation of the infinite case.
    # NUM_SIMULATIONS: Number of trials to run. More trials lead to a more
    #                  stable and accurate average.
    N = 2000
    NUM_SIMULATIONS = 50

    total_stars_simulated = 0
    total_constellations_found = 0

    print(f"Running {NUM_SIMULATIONS} simulations with {N} stars each...")

    for i in range(NUM_SIMULATIONS):
        # 1. Generate N random points in a [0,1] x [0,1] square
        points = np.random.rand(N, 2)

        # 2. For each point, find its nearest neighbor using toroidal distance
        # This creates the directed edges of the graph
        adj = [set() for _ in range(N)]
        for j in range(N):
            min_dist_sq = float('inf')
            nn_index = -1
            p1 = points[j]
            for k in range(N):
                if j == k:
                    continue
                
                p2 = points[k]
                # Calculate squared toroidal distance (faster than sqrt)
                delta = np.abs(p1 - p2)
                # If distance is > 0.5, wrap around the torus
                delta = np.where(delta > 0.5, 1.0 - delta, delta)
                dist_sq = np.sum(delta**2)

                if dist_sq < min_dist_sq:
                    min_dist_sq = dist_sq
                    nn_index = k
            
            # 3. Build the undirected graph adjacency list
            # An edge from j to its nearest neighbor nn_index means they are
            # in the same component.
            adj[j].add(nn_index)
            adj[nn_index].add(j)

        # For component counting, convert sets to lists
        undirected_adj = [list(s) for s in adj]

        # 4. Count connected components using iterative DFS
        visited = [False] * N
        num_components = 0
        for j in range(N):
            if not visited[j]:
                num_components += 1
                stack = [j]
                visited[j] = True
                while stack:
                    u = stack.pop()
                    for v in undirected_adj[u]:
                        if not visited[v]:
                            visited[v] = True
                            stack.append(v)
        
        total_stars_simulated += N
        total_constellations_found += num_components
        
        # Optional: print progress
        # print(f"Simulation {i+1}/{NUM_SIMULATIONS} done. Found {num_components} constellations.")


    # 5. Calculate the final average size
    average_size = total_stars_simulated / total_constellations_found

    print("\n--- Simulation Results ---")
    print(f"Total stars simulated: {total_stars_simulated}")
    print(f"Total constellations found: {total_constellations_found}")
    print(f"Average number of stars per constellation = {total_stars_simulated} / {total_constellations_found} = {average_size:.4f}")
    
    # The final answer is the calculated average size.
    return average_size

# Run the simulation and print the final answer in the required format.
final_answer = solve_constellation_problem()
# The <<<...>>> format is for automated result extraction.
print(f"\n<<< {final_answer:.4f} >>>")
