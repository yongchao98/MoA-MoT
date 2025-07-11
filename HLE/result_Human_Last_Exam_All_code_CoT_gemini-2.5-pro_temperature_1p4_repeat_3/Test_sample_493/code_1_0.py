import numpy as np
from scipy.spatial import cKDTree

def calculate_average_constellation_size(num_stars, box_size=1.0):
    """
    Simulates stars in a 2D patch of sky and calculates the average
    number of stars per constellation.

    Args:
        num_stars (int): The total number of stars to simulate.
        box_size (float): The side length of the simulation area.

    Returns:
        None. Prints the result of the calculation.
    """
    if num_stars < 2:
        print("Cannot form constellations with fewer than 2 stars.")
        return

    # 1. Generate N stars (points) in a square with periodic boundary conditions.
    # We use a uniform random distribution.
    points = np.random.rand(num_stars, 2) * box_size

    # 2. Find the nearest neighbor for each star using cKDTree.
    # The `boxsize` argument tells the cKDTree to use periodic (toroidal) boundary conditions.
    tree = cKDTree(points, boxsize=[box_size, box_size])
    
    # We query for the 2 nearest neighbors because the first (at index 0) is the point itself.
    # The actual nearest neighbor is at index 1.
    distances, indices = tree.query(points, k=2)

    # nearest_neighbors[i] will store the index of the nearest neighbor to star i.
    nearest_neighbors = indices[:, 1]

    # 3. Count the number of mutual nearest-neighbor (MNN) pairs.
    # This count is equal to the number of constellations.
    num_constellations = 0
    
    # Iterate through each star to check if it's part of an MNN pair.
    for i in range(num_stars):
        j = nearest_neighbors[i]
        # To avoid double-counting pairs (i,j) and (j,i), we only check pairs where i < j.
        if i < j:
            # Check for mutuality: is the nearest neighbor of star j also star i?
            if nearest_neighbors[j] == i:
                num_constellations += 1

    # 4. Calculate and print the average number of stars per constellation.
    if num_constellations > 0:
        avg_size = num_stars / num_constellations
        print(f"Total stars simulated (N): {num_stars}")
        print(f"Number of constellations found (MNN pairs): {num_constellations}")
        print(f"Average stars per constellation (N / MNN pairs): {num_stars} / {num_constellations} = {avg_size}")
    else:
        print(f"No mutual nearest-neighbor pairs found among {num_stars} stars.")


if __name__ == '__main__':
    # A large number of stars provides a more accurate statistical average.
    N_STARS = 100000
    calculate_average_constellation_size(N_STARS)
    
# The theoretical value from stochastic geometry is approximately 3.218
<<<3.218>>>