import numpy as np

def generate_stretched_grid(num_points, stretching_factor):
    """
    Demonstrates grid stretching using the tanh function.

    Args:
        num_points (int): The number of grid points to generate.
        stretching_factor (float): A factor to control the intensity of the stretching.
                                   Higher values lead to more clustering at the ends.
    """
    # 1. Create a uniform grid in "computational space" from -1 to 1.
    uniform_grid = np.linspace(-1.0, 1.0, num_points)

    # 2. Apply the tanh stretching function.
    # The division by np.tanh(stretching_factor) normalizes the output to remain in the [-1, 1] range.
    stretched_grid = np.tanh(stretching_factor * uniform_grid) / np.tanh(stretching_factor)

    # 3. Print the results for comparison.
    print("--- Grid Stretching Demonstration using tanh ---")
    print(f"\nOriginal Uniform Grid ({num_points} points):")
    # np.round is used for cleaner output
    print(np.round(uniform_grid, 4))

    print(f"\nTanh-Stretched Grid (Stretching Factor: {stretching_factor}):")
    print("Notice how the points are now clustered towards -1.0 and 1.0.")
    print(np.round(stretched_grid, 4))

# --- Parameters ---
# Number of points in our 1D grid
number_of_points = 11
# How much to cluster points at the ends. Try changing this value (e.g., to 1.0 or 3.0)
stretch_k = 2.0

generate_stretched_grid(number_of_points, stretch_k)