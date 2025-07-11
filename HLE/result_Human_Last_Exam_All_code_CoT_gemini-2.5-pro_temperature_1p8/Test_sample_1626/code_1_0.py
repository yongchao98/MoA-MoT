import numpy as np

def demonstrate_tanh_grid_clustering():
    """
    Demonstrates how the tanh function is used to cluster grid points.
    """
    # Number of grid points
    num_points = 11

    # Strength of the clustering. Higher values mean more clustering.
    # A value of 0 would result in no clustering (uniform grid).
    clustering_strength = 2.5

    # 1. Create a uniform grid in computational space from -1 to 1.
    # This represents our initial, evenly spaced grid.
    uniform_grid = np.linspace(-1, 1, num_points)

    # 2. Apply the tanh transformation to get the clustered grid in physical space.
    # This formula maps [-1, 1] to [-1, 1] while clustering points around 0.
    # We use np.tanh(clustering_strength) in the denominator to normalize the output
    # so that the boundaries remain at -1 and 1.
    clustered_grid = np.tanh(clustering_strength * uniform_grid) / np.tanh(clustering_strength)

    print("Original Uniform Grid Points:")
    # We print each number from the original grid
    for point in uniform_grid:
        print(f"{point:8.4f}", end=" ")
    print("\n")

    print(f"Grid Points After Tanh Transformation (Clustering Strength = {clustering_strength}):")
    # We print each number from the transformed grid, which results from the function's application
    for point in clustered_grid:
        print(f"{point:8.4f}", end=" ")
    print("\n")
    print("Notice how the points in the second grid are clustered more densely around 0.0000.")


demonstrate_tanh_grid_clustering()