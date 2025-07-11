import numpy as np

def demonstrate_tanh_grid_spacing():
    """
    Demonstrates how the hyperbolic tangent (tanh) function is used to
    govern grid spacing and create a clustered grid from a uniform one.
    """
    print("--- Demonstrating Grid Spacing with the tanh Function ---")

    # 1. Define parameters
    num_points = 21  # Total number of grid points
    # The clustering parameter 'k' controls the intensity of clustering.
    # A value close to 0 results in a nearly uniform grid.
    # A large value results in strong clustering towards the center.
    k_clustering = 2.5

    # 2. Create a uniform grid
    # This represents the initial, evenly-spaced points in a computational space, e.g., from -1 to 1.
    uniform_grid = np.linspace(-1, 1, num_points)

    # 3. Apply the tanh mapping function
    # This is a standard formula for creating a stretched grid.
    # The tanh function smoothly remaps the points from the uniform grid.
    # Dividing by tanh(k) normalizes the output to ensure the mapped grid still spans from -1 to 1.
    clustered_grid = np.tanh(k_clustering * uniform_grid) / np.tanh(k_clustering)

    # 4. Print the results for comparison
    print(f"\nComparing a uniform grid with a tanh-clustered grid for {num_points} points.")
    print(f"The clustering factor (k) is: {k_clustering}\n")

    print("Original Uniform Grid Points:")
    print(np.round(uniform_grid, 4))
    print("\nResulting Clustered Grid Points (using tanh):")
    print(np.round(clustered_grid, 4))

    print("\nAnalysis:")
    print("Notice how the points in the 'Clustered Grid' are more densely packed near the center (0.0)")
    print("and become progressively more spread out towards the edges (-1.0 and 1.0).")
    print("This smooth control over grid spacing is why 'tanh' is the key function for this purpose.")


# Execute the demonstration
demonstrate_tanh_grid_spacing()