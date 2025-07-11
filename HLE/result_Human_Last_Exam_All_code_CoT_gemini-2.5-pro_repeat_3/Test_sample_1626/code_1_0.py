import numpy as np

def demonstrate_tanh_spacing():
    """
    Demonstrates how the tanh function is used for grid spacing.
    It creates a uniform grid and then a tanh-stretched grid to show
    how points are clustered towards the boundaries.
    """
    num_points = 11  # Number of grid points
    stretching_factor = 2.0  # A factor > 0. Higher values mean more stretching.

    # 1. Create a uniform grid in the computational domain [-1, 1]
    uniform_grid = np.linspace(-1, 1, num_points)

    # 2. Apply the tanh transformation to create a stretched physical grid
    # The equation used is: x_physical = tanh(stretching_factor * x_uniform) / tanh(stretching_factor)
    # This keeps the grid within the [-1, 1] domain but clusters points near -1 and 1.
    denominator = np.tanh(stretching_factor)
    stretched_grid = np.tanh(stretching_factor * uniform_grid) / denominator

    print("Demonstration of Tanh Grid Spacing")
    print("-" * 40)
    print(f"Number of points: {num_points}")
    print(f"Stretching factor: {stretching_factor}")
    print("-" * 40)

    print("Uniform Grid Points:")
    # Using a loop to print each number clearly
    for point in uniform_grid:
        print(f"{point:8.4f}", end=" ")
    print("\n")

    print("Tanh-Stretched Grid Points (clustered at ends):")
    # This is the final result of the equation, with each number outputted
    for point in stretched_grid:
        print(f"{point:8.4f}", end=" ")
    print("\n")
    print("-" * 40)
    print("Notice how the points in the stretched grid are closer together near -1.0 and 1.0.")

demonstrate_tanh_spacing()