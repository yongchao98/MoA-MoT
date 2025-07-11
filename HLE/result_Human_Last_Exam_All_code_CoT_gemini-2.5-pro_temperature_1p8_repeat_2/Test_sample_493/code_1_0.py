import numpy as np
from scipy.spatial import KDTree

def calculate_average_constellation_size():
    """
    Simulates a nearest-neighbor graph for a uniform distribution of stars
    to find the average size of a connected component (constellation).
    """
    # Number of stars to simulate. A larger number gives a more accurate estimate.
    N = 20000

    # Step 1: Generate N random points (stars) in a [0,1]x[0,1] square.
    # This represents a uniform distribution in the sky.
    points = np.random.rand(N, 2)

    # Step 2: Handle toroidal (wrap-around) boundaries to avoid edge effects.
    # We do this by tiling the points in a 3x3 grid around the central unit square.
    tiled_points = []
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            tiled_points.append(points + np.array([dx, dy]))
    tiled_points = np.vstack(tiled_points)

    # Step 3: Use a KD-tree for efficient nearest-neighbor search.
    # The KD-tree is built on the tiled set of points.
    kdtree = KDTree(tiled_points)

    # Step 4: For each original point, find its nearest neighbor from the tiled set.
    # We query for k=2 because the closest point (at index 0) to a point `p` is `p` itself.
    # The actual nearest neighbor is the second closest (at index 1).
    # We use `workers=-1` to utilize all available CPU cores for speed.
    distances, indices = kdtree.query(points, k=2, workers=-1)

    # The index of the nearest neighbor for point `i` is `indices[i, 1]`.
    # This index is for the tiled array, so we map it back to the original 0 to N-1 range.
    nearest_neighbor_idx = indices[:, 1] % N

    # Step 5: Use a Disjoint Set Union (DSU) data structure to find connected components.
    # This is an efficient way to group vertices connected by edges.
    parent = list(range(N))
    def find_set(v):
        if v == parent[v]:
            return v
        # Path compression for optimization
        parent[v] = find_set(parent[v])
        return parent[v]

    def unite_sets(a, b):
        a = find_set(a)
        b = find_set(b)
        if a != b:
            # Union by rank/size could be used, but this is sufficient.
            parent[b] = a

    # For each star `i` and its nearest neighbor `j`, unite their sets.
    # Since "constellations" are connected subgraphs, we ignore edge direction
    # and unite the sets for (star, its_neighbor).
    for i in range(N):
        j = nearest_neighbor_idx[i]
        unite_sets(i, j)

    # Step 6: Count the number of distinct constellations.
    # This is the number of unique roots in our DSU structure.
    num_components = len(set(find_set(i) for i in range(N)))

    # Step 7: Calculate the average constellation size.
    average_size = N / num_components

    # Print the results.
    print(f"Simulation based on a model with {N} stars.")
    print("A constellation is a connected component in a nearest-neighbor graph.")
    print("The final equation is: Average Size = Total Stars / Number of Constellations")
    print(f"Average Size = {N} / {num_components}")
    print(f"The average number of stars per constellation is: {average_size:.4f}")

# Run the simulation and print the output.
calculate_average_constellation_size()