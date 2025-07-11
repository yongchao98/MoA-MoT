import numpy as np
from scipy.spatial import KDTree

def find_average_constellation_size():
    """
    Simulates the constellation model and calculates the average size.

    This function performs the following steps:
    1. Generates N random stars in a 2D plane.
    2. For each star, finds its nearest neighbor efficiently using a KDTree.
    3. Builds the underlying undirected graph where stars are nodes and a
       star is connected to its nearest neighbor.
    4. Counts the number of connected components (constellations) in this graph.
    5. Calculates the average number of stars per constellation.
    """
    # Step 1: Model the Stars
    # A large number of stars gives a better approximation of the theoretical value.
    num_stars = 20000
    
    # Use a fixed seed for reproducibility of the random star placement.
    np.random.seed(42)
    stars = np.random.rand(num_stars, 2)

    # Step 2: Identify Nearest Neighbors
    # We use a KDTree for efficient nearest neighbor searching (O(N log N)).
    # We query for the 2 nearest neighbors because the first one (at index 0)
    # is always the point itself.
    kdtree = KDTree(stars)
    _, indices = kdtree.query(stars, k=2)
    
    # The nearest neighbor of star `i` is at `indices[i, 1]`.
    nearest_neighbors = indices[:, 1]

    # Step 3: Form Constellations (build undirected graph)
    # We create an adjacency list for the undirected graph.
    # If star `i`'s nearest neighbor is `j`, we connect both i-j and j-i.
    adj = [[] for _ in range(num_stars)]
    for i in range(num_stars):
        j = nearest_neighbors[i]
        adj[i].append(j)
        adj[j].append(i)

    # Step 4: Count Constellations (connected components)
    visited = [False] * num_stars
    num_components = 0
    for i in range(num_stars):
        if not visited[i]:
            num_components += 1
            # Start a traversal (DFS) from node i to find all nodes in its component.
            stack = [i]
            visited[i] = True
            while stack:
                u = stack.pop()
                # Use a set for the adjacency list to handle potential duplicates efficiently
                for v in set(adj[u]):
                    if not visited[v]:
                        visited[v] = True
                        stack.append(v)

    # Step 5: Calculate and print the Average
    if num_components > 0:
        average_size = num_stars / num_components
        # The final output prints the equation as requested
        print(f"Total stars simulated: {num_stars}")
        print(f"Constellations found: {num_components}")
        print(f"Average stars per constellation calculation:")
        print(f"{num_stars} / {num_components} = {average_size}")
    else:
        print("No components found.")

if __name__ == '__main__':
    find_average_constellation_size()