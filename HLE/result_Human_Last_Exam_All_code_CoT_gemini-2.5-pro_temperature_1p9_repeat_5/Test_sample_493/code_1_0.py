import numpy as np
from scipy.spatial import cKDTree

def calculate_average_constellation_size(num_stars=20000):
    """
    Simulates the constellation problem to find the average number of stars
    per connected component.

    Args:
        num_stars (int): The number of stars to include in the simulation.
                         A larger number provides a more accurate estimate.

    Returns:
        None. Prints the result of the simulation.
    """
    if num_stars < 2:
        print("Need at least 2 stars for the simulation.")
        return

    # Step 1: Generate N random points (stars) in a 2D plane (unit square)
    points = np.random.rand(num_stars, 2)

    # Step 2: Find the nearest neighbor for each point using a k-d tree.
    # We query for k=2 because the closest point to any point is the point itself.
    # The second closest is its nearest neighbor.
    kdtree = cKDTree(points)
    distances, indices = kdtree.query(points, k=2)
    nearest_neighbor_indices = indices[:, 1]

    # Step 3: Build an undirected graph using an adjacency list.
    # An edge exists between A and B if B is the nearest neighbor of A (or vice versa).
    adj_list = {i: [] for i in range(num_stars)}
    for i in range(num_stars):
        neighbor = nearest_neighbor_indices[i]
        # Add an undirected edge between the point and its nearest neighbor
        adj_list[i].append(neighbor)
        adj_list[neighbor].append(i)

    # Step 4: Count the number of connected components (constellations) using DFS.
    visited = [False] * num_stars
    num_components = 0
    for i in range(num_stars):
        # If we find a star that hasn't been visited, it's part of a new component.
        if not visited[i]:
            num_components += 1
            # Start a Depth-First Search (DFS) to find all stars in this component.
            stack = [i]
            visited[i] = True
            while stack:
                current_star = stack.pop()
                for connected_star in adj_list[current_star]:
                    if not visited[connected_star]:
                        visited[connected_star] = True
                        stack.append(connected_star)

    # Step 5: Calculate the average and print the result.
    if num_components > 0:
        average_size = num_stars / num_components
        print(f"Simulation based on {num_stars} stars.")
        print("Total stars / Total constellations = Average stars per constellation")
        print(f"{num_stars} / {num_components} = {average_size:.4f}")
    else:
        print("No constellations found.")

if __name__ == '__main__':
    # Running the simulation with 20,000 stars for a good estimate.
    # The result may vary slightly on each execution due to randomness.
    calculate_average_constellation_size()