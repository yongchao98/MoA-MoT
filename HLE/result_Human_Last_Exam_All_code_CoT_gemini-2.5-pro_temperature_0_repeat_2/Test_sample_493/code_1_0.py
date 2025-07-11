import numpy as np
from scipy.spatial import cKDTree
from collections import deque

def solve_constellation_problem():
    """
    Simulates the constellation model to find the average number of stars per constellation.
    """
    # Step 1: Generate a large number of stars (points) uniformly in a 2D space.
    # A larger N gives a more accurate approximation of the infinite case.
    N = 50000
    
    # Generate N points in a 2D unit square. Boundary effects are minimized for large N.
    points = np.random.rand(N, 2)

    # Step 2: For each star, find its nearest neighbor to build the graph.
    # We use a k-d tree for efficient nearest neighbor searching (O(N log N)).
    # We query for k=2 because the nearest point to any point is itself.
    tree = cKDTree(points)
    # The second column of the result gives the index of the nearest neighbor.
    distances, indices = tree.query(points, k=2)
    nearest_neighbors = indices[:, 1]

    # Step 3: Identify constellations (connected components) in the undirected graph.
    # We build an adjacency list for the undirected version of the graph.
    # If i's nearest neighbor is j, we consider them connected.
    undirected_graph = [[] for _ in range(N)]
    for i in range(N):
        j = nearest_neighbors[i]
        undirected_graph[i].append(j)
        undirected_graph[j].append(i)

    # Use Breadth-First Search (BFS) to count the number of connected components.
    visited = [False] * N
    num_components = 0
    for i in range(N):
        if not visited[i]:
            # Found a new, unvisited component.
            num_components += 1
            q = deque([i])
            visited[i] = True
            # Traverse all nodes in this component to mark them as visited.
            while q:
                u = q.popleft()
                for v in undirected_graph[u]:
                    if not visited[v]:
                        visited[v] = True
                        q.append(v)

    # Step 4: Calculate the average size and print the result.
    if num_components > 0:
        average_size = N / num_components
        print("This simulation models stars as uniformly distributed points in a 2D plane.")
        print("A constellation is a group of stars connected by nearest-neighbor links.")
        print("\nCalculating the average number of stars per constellation:")
        print(f"Total number of stars simulated: {N}")
        print(f"Number of constellations found: {num_components}")
        print("\nFinal Equation:")
        print(f"{N} / {num_components} = {average_size:.4f}")
    else:
        print("No components found.")

# Run the simulation.
# Note: The result will vary slightly each time due to the random placement of stars.
solve_constellation_problem()