import numpy as np
from sklearn.neighbors import NearestNeighbors
from collections import deque

def find_average_constellation_size():
    """
    This function simulates the star constellation model to find the average
    number of stars per constellation.
    """
    # Step 1: Define simulation parameters. A large number of stars provides a
    # more accurate approximation of the theoretical average.
    num_stars = 20000
    
    # Use a fixed seed for the random number generator to ensure the result is reproducible.
    np.random.seed(42)

    # Step 2: Generate random star positions uniformly in a 2D unit square.
    stars = np.random.rand(num_stars, 2)

    # Step 3: Find the single nearest neighbor for each star.
    # We ask for 2 neighbors because the closest point to any given star is the star itself.
    # The second closest is its nearest distinct neighbor.
    nbrs = NearestNeighbors(n_neighbors=2, algorithm='auto').fit(stars)
    distances, indices = nbrs.kneighbors(stars)
    
    # The nearest neighbor for star `i` is located at `indices[i, 1]`.
    nearest_neighbor_indices = indices[:, 1]

    # Step 4: Build the undirected graph. A constellation is a connected component.
    # We use an adjacency list to represent the graph.
    adj_list = [[] for _ in range(num_stars)]
    for i in range(num_stars):
        j = nearest_neighbor_indices[i]
        # Add an edge in both directions to make the graph undirected.
        adj_list[i].append(j)
        adj_list[j].append(i)

    # Step 5: Count the number of connected components (constellations) using BFS.
    visited = [False] * num_stars
    num_constellations = 0
    for i in range(num_stars):
        if not visited[i]:
            num_constellations += 1
            # Start a traversal (BFS) to find all stars in this component.
            queue = deque([i])
            visited[i] = True
            while queue:
                current_star = queue.popleft()
                for neighbor in adj_list[current_star]:
                    if not visited[neighbor]:
                        visited[neighbor] = True
                        queue.append(neighbor)
    
    # Step 6: Calculate the average size and print the result.
    average_size = num_stars / num_constellations
    
    print(f"Total number of stars simulated: {num_stars}")
    print(f"Total number of constellations found: {num_constellations}")
    print("Average number of stars per constellation:")
    print(f"{num_stars} / {num_constellations} = {average_size}")


if __name__ == '__main__':
    find_average_constellation_size()