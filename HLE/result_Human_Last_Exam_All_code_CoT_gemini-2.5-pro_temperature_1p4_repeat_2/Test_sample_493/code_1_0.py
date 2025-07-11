import numpy as np
from scipy.spatial import KDTree

def find_average_constellation_size(num_stars, random_seed=42):
    """
    Simulates the formation of constellations and calculates their average size.

    Args:
        num_stars (int): The number of stars to include in the simulation.
        random_seed (int): A seed for the random number generator for reproducibility.
    """
    # Set a seed for reproducibility of the random star positions
    np.random.seed(random_seed)

    # 1. Model the Stars: Generate N random stars in a 2D unit square.
    # Each row in the 'stars' array is a star with [x, y] coordinates.
    stars = np.random.rand(num_stars, 2)

    # 2. Find Nearest Neighbors: Use a k-d tree for efficient lookup.
    # We query for the 2 nearest neighbors because the closest neighbor to any
    # point is the point itself (at distance 0). The actual nearest neighbor
    # is the second one in the returned list.
    kdtree = KDTree(stars)
    # indices[i, 1] will hold the index of the nearest neighbor to star i.
    _, indices = kdtree.query(stars, k=2)
    nearest_neighbors = indices[:, 1]

    # 3. Identify Constellations: Find connected components.
    # We build an undirected graph adjacency list to find the components.
    # If star 'i' is connected to 'j', we add an edge in both directions.
    adj = [[] for _ in range(num_stars)]
    for i in range(num_stars):
        j = nearest_neighbors[i]
        adj[i].append(j)
        adj[j].append(i)

    # Use a Depth-First Search (DFS) based approach to count components.
    visited = [False] * num_stars
    num_components = 0
    for i in range(num_stars):
        if not visited[i]:
            num_components += 1
            # Start a traversal from this unvisited star
            stack = [i]
            visited[i] = True
            while stack:
                u = stack.pop()
                for v in adj[u]:
                    if not visited[v]:
                        visited[v] = True
                        stack.append(v)

    # 4. Calculate the Average
    if num_components > 0:
        average_size = num_stars / num_components
    else:
        average_size = 0

    # Output the numbers used in the final equation, as requested.
    print(f"Simulation run with {num_stars} stars.")
    print(f"Found {num_components} separate constellations.")
    print(f"Average stars per constellation = Total Stars / Number of Constellations")
    print(f"Final Calculation: {num_stars} / {num_components} = {average_size:.4f}")

if __name__ == '__main__':
    # Using a large number of stars gives a more accurate and stable result.
    NUMBER_OF_STARS = 50000
    find_average_constellation_size(NUMBER_OF_STARS)