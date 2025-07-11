import numpy as np
from scipy.spatial import cKDTree

def solve_constellations():
    """
    Simulates the constellation model to find the average number of stars
    per constellation.
    """
    # --- Step 1: Model the Stars ---
    # N should be large enough to get a stable statistical average.
    N_points = 20000
    # Generate points in a 1x1 square for our main viewing area.
    points = np.random.rand(N_points, 2)

    # --- Step 2: Handle Boundaries (Periodic/Tiling Method) ---
    # Create a 3x3 grid of the points to simulate a toroidal space.
    # This ensures that points near an edge can find neighbors on the "other side".
    tiled_points = []
    original_indices = []
    for i in range(-1, 2):
        for j in range(-1, 2):
            offset = np.array([i, j])
            tiled_points.append(points + offset)
            # We need to remember which original point each tiled point corresponds to.
            original_indices.append(np.arange(N_points))

    all_points_tiled = np.vstack(tiled_points)
    all_original_indices = np.hstack(original_indices)

    # --- Step 3: Find Nearest Neighbors ---
    # Build a k-d tree on the tiled set of points for fast searching.
    # A cKDTree is a C-implemented version, which is faster.
    kdtree = cKDTree(all_points_tiled)

    # For each original point, find its nearest neighbor in the tiled universe.
    # We query for k=2 because the point itself is the nearest (distance 0).
    # Its true nearest neighbor will be the second one returned.
    distances, indices = kdtree.query(points, k=2)
    
    # Get the indices of the nearest neighbors. These are indices into 'all_points_tiled'.
    nn_tiled_indices = indices[:, 1]
    # Map these back to the original point indices (0 to N-1).
    nn_original_indices = all_original_indices[nn_tiled_indices]

    # --- Step 4 & 5: Form and Count Constellations (Union-Find) ---
    # A simple implementation of a Union-Find data structure.
    parent = list(range(N_points))
    num_components = N_points

    def find_set(v):
        if v == parent[v]:
            return v
        parent[v] = find_set(parent[v]) # Path compression
        return parent[v]

    def unite_sets(a, b):
        nonlocal num_components
        a = find_set(a)
        b = find_set(b)
        if a != b:
            parent[b] = a
            num_components -= 1
    
    # For each star, create a link to its nearest neighbor, forming the components.
    for i in range(N_points):
        unite_sets(i, nn_original_indices[i])

    # --- Step 6: Calculate the Average ---
    # The average size is the total number of stars divided by the number of components.
    average_size = N_points / num_components

    # --- Final Output ---
    print(f"Simulation run with N = {N_points} stars.")
    print(f"Found C = {num_components} distinct constellations.")
    print(f"The average number of stars per constellation (N / C) is:")
    # The final equation asked for in the prompt.
    print(f"{N_points} / {num_components} = {average_size}")
    return average_size

# Run the simulation and capture the final answer.
final_answer = solve_constellations()
# Use <<<>>> format for the final numerical answer.
print(f"\n<<<__{final_answer:.4f}__>>>")
