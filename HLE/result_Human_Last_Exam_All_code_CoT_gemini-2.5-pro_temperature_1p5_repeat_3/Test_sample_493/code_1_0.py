import numpy as np
from scipy.spatial import cKDTree

def solve_constellation_problem():
    """
    Simulates the constellation problem to find the average number of stars per constellation.
    """
    # Number of stars in the simulation. A larger number gives a more accurate result.
    num_stars = 20000

    # Step 1: Generate star positions
    # We generate N points uniformly distributed in a 2D square.
    print(f"Simulating a sky with {num_stars} stars...")
    points = np.random.rand(num_stars, 2)

    # Step 2: Find the nearest neighbor for each star
    # We use a k-d tree for efficient searching. We query for k=2 because
    # the closest point to any point (at distance 0) is the point itself.
    tree = cKDTree(points)
    # indices[:, 0] will be the point itself, indices[:, 1] will be its nearest neighbor.
    _, indices = tree.query(points, k=2)
    nearest_neighbors = indices[:, 1]

    # Step 3: Determine constellations using a Union-Find data structure.
    # Initially, each star is in its own set (constellation).
    parent = list(range(num_stars))
    
    def find_set(v):
        """Find the representative of the set containing element v with path compression."""
        if v == parent[v]:
            return v
        parent[v] = find_set(parent[v])
        return parent[v]

    def unite_sets(a, b):
        """Unite the sets containing elements a and b."""
        a_root = find_set(a)
        b_root = find_set(b)
        if a_root != b_root:
            parent[b_root] = a_root
            return True # A merge happened
        return False # Already in the same set

    # For each star i, it's connected to its nearest neighbor j.
    # We unite their sets.
    for i in range(num_stars):
        j = nearest_neighbors[i]
        unite_sets(i, j)

    # Step 4: Count the number of unique constellations.
    # This is the number of root nodes in our Union-Find structure.
    num_constellations = 0
    for i in range(num_stars):
        if parent[i] == i:
            num_constellations += 1

    # Step 5: Calculate and print the final result.
    if num_constellations > 0:
        average_size = num_stars / num_constellations
        print(f"\nCalculation:")
        print(f"Total Stars / Number of Constellations = Average Size")
        print(f"{num_stars} / {num_constellations} = {average_size:.4f}")
    else:
        print("No constellations found.")

solve_constellation_problem()