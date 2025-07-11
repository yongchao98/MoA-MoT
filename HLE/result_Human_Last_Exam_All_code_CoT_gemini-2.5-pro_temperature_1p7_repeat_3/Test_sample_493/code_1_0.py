import numpy as np
import sys

def solve_constellation_problem():
    """
    Simulates star constellations to find the average number of stars per constellation.
    
    This function performs a single, large-scale simulation.
    1. Generates N random points (stars).
    2. For each star, finds its nearest neighbor.
    3. Constructs an undirected graph from these connections.
    4. Counts the number of connected components (constellations).
    5. Calculates and prints the average constellation size.
    """
    # Number of stars in the simulation. A larger number gives a more accurate result.
    num_stars = 5000
    
    # Increase recursion limit for DFS on potentially long, chain-like components.
    sys.setrecursionlimit(num_stars + 50)

    # 1. Generate random star positions in a 2D unit square.
    points = np.random.rand(num_stars, 2)

    # 2. For each star, find its nearest neighbor.
    #    This creates a directed graph where successors[i] is the neighbor of i.
    successors = np.zeros(num_stars, dtype=int)
    for i in range(num_stars):
        # Calculate squared Euclidean distances from point i to all other points.
        # Using squared distance is faster as it avoids the sqrt operation.
        diff = points - points[i]
        dist_sq = np.sum(diff**2, axis=1)
        
        # Set distance to self to infinity to avoid picking the point itself.
        dist_sq[i] = np.inf
        
        # Find the index of the minimum distance.
        successors[i] = np.argmin(dist_sq)

    # 3. Build the undirected graph from the directed nearest neighbor links.
    #    An edge (u, v) exists if v is the nearest neighbor of u.
    #    We use an adjacency list to represent the graph.
    adj = [[] for _ in range(num_stars)]
    for i in range(num_stars):
        j = successors[i]
        # Add edge in both directions for the undirected graph.
        adj[i].append(j)
        adj[j].append(i)

    # 4. Count connected components (constellations) using Depth First Search (DFS).
    visited = [False] * num_stars
    component_count = 0

    def dfs(u):
        """Recursive DFS to traverse a component."""
        visited[u] = True
        for v in adj[u]:
            if not visited[v]:
                dfs(v)

    for i in range(num_stars):
        if not visited[i]:
            component_count += 1
            dfs(i)

    # 5. Calculate the average constellation size.
    avg_size = num_stars / component_count if component_count > 0 else 0
    
    # Final Output
    print("This simulation models stars as points and connects each to its nearest neighbor.")
    print("The average size of a constellation is the total number of stars divided by the number of connected components.\n")
    print("Equation from this simulation run:")
    print(f"Total Stars / Number of Constellations = Average Stars per Constellation")
    print(f"{num_stars} / {component_count} = {avg_size:.4f}")

# Run the simulation and print the results.
solve_constellation_problem()