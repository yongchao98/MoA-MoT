import numpy as np
from scipy.spatial import cKDTree

def solve_constellation_problem():
    """
    Simulates the constellation model to find the average number of stars
    per constellation.
    """
    # Step 1: Generate a large number of stars on a 2D torus for simulation.
    # Using a fixed seed for reproducibility of the random process.
    np.random.seed(42)
    num_stars = 20000
    # The points are generated in a unit square [0,1]x[0,1]
    points = np.random.rand(num_stars, 2)

    # Step 2: Find the nearest neighbor for each star.
    # We use a cKDTree with periodic boundary conditions (boxsize=[1.0, 1.0])
    # to simulate a torus and avoid edge effects.
    # We query for k=2 because the closest point to any point is itself.
    tree = cKDTree(points, boxsize=[1.0, 1.0])
    _, indices = tree.query(points, k=2)
    # The nearest neighbor of each point is in the second column.
    nearest_neighbors = indices[:, 1]

    # Step 3: Build the adjacency list for the UNDIRECTED graph.
    # An edge exists between i and j if i is NN of j OR j is NN of i.
    adj = [[] for _ in range(num_stars)]
    for i in range(num_stars):
        j = nearest_neighbors[i]
        adj[i].append(j)
        adj[j].append(i)

    # Step 4: Count the connected components (constellations) using BFS.
    visited = [False] * num_stars
    num_components = 0
    for i in range(num_stars):
        if not visited[i]:
            num_components += 1
            # Start a Breadth-First Search (BFS) to find all stars in this component
            q = [i]
            visited[i] = True
            head = 0
            while head < len(q):
                u = q[head]
                head += 1
                for v in adj[u]:
                    if not visited[v]:
                        visited[v] = True
                        q.append(v)
    
    # Step 5: Calculate and print the average constellation size.
    # The average size is the total number of stars divided by the number of components.
    average_size = num_stars / num_components
    
    print("--- Simulation Results ---")
    print(f"Total number of stars (N): {num_stars}")
    print(f"Number of constellations found (C): {num_components}")
    print("\n--- Final Calculation ---")
    print(f"Average stars per constellation = N / C")
    print(f"= {num_stars} / {num_components}")
    print(f"= {average_size}")

solve_constellation_problem()