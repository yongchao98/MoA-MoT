import numpy as np
from scipy.spatial import cKDTree

def solve_constellation_problem():
    """
    Simulates the star constellation model to find the average number of stars
    per constellation.
    """
    # Number of stars to simulate. A larger number gives a more accurate result.
    num_stars = 100000
    
    # The side length of the 2D patch of sky. We can set this to 1.0.
    box_size = 1.0

    # 1. Generate random star positions in a 2D square.
    # np.random.rand(num_stars, 2) creates an array of shape (num_stars, 2)
    # with random coordinates between 0.0 and 1.0.
    stars = np.random.rand(num_stars, 2) * box_size

    # 2. Find the nearest neighbor for each star using a k-d tree.
    # We use cKDTree because it's highly optimized for this task and can
    # handle periodic (toroidal) boundary conditions via the `boxsize` parameter.
    tree = cKDTree(stars, boxsize=[box_size, box_size])
    
    # For each star, we query for its 2 nearest neighbors. The first (k=1) will
    # be the star itself, so the second (k=2) is its actual nearest neighbor.
    _, nn_indices = tree.query(stars, k=2)
    
    # The array of nearest neighbors for each star i is nn_indices[i, 1].
    nearest_neighbors = nn_indices[:, 1]

    # 3. Use a Disjoint Set Union (DSU) data structure to count connected components.
    # Each component represents one constellation.
    parent = list(range(num_stars))
    num_constellations = num_stars

    def find_set(v):
        # Finds the representative (root) of the set containing element v
        if v == parent[v]:
            return v
        # Path compression for efficiency
        parent[v] = find_set(parent[v])
        return parent[v]

    def unite_sets(a, b):
        # Merges the sets containing elements a and b
        nonlocal num_constellations
        a_root = find_set(a)
        b_root = find_set(b)
        if a_root != b_root:
            parent[b_root] = a_root
            num_constellations -= 1

    # For each star, we form a connection to its nearest neighbor,
    # merging their sets if they are not already connected.
    for i in range(num_stars):
        unite_sets(i, nearest_neighbors[i])

    # 4. Calculate the average size of a constellation.
    avg_size = num_stars / num_constellations
    
    # Output the components of the final calculation as requested.
    print(f"Total number of stars simulated: {num_stars}")
    print(f"Number of constellations found: {num_constellations}")
    print(f"Average stars per constellation: {num_stars} / {num_constellations} = {avg_size}")
    
    # The final answer in the required format.
    print(f"\n<<<{avg_size}>>>")

# Run the simulation and print the results
solve_constellation_problem()