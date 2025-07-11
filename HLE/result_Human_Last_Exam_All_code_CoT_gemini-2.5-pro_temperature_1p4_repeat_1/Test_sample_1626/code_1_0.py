import numpy as np

def generate_stretched_grid(num_points, stretching_factor):
    """
    Generates a grid stretched using the tanh function.

    Args:
        num_points (int): The number of grid points.
        stretching_factor (float): Controls the intensity of the clustering.
                                   A value of 0 means no stretching.
                                   Higher values (e.g., 2, 3) lead to more intense clustering at the ends.
    """
    if stretching_factor < 0:
        raise ValueError("Stretching factor must be non-negative.")

    # 1. Create a uniform grid in computational space from -1 to 1
    uniform_grid = np.linspace(-1.0, 1.0, num_points)

    # 2. Apply the tanh stretching function
    # The division by np.tanh(stretching_factor) normalizes the output
    # to ensure the domain remains [-1, 1].
    # A stretching_factor of 0 would cause division by zero, so we handle it.
    if stretching_factor == 0:
        stretched_grid = uniform_grid
    else:
        stretched_grid = np.tanh(stretching_factor * uniform_grid) / np.tanh(stretching_factor)

    print("--- Grid Generation using tanh Function ---")
    print(f"Number of points: {num_points}")
    print(f"Stretching factor: {stretching_factor}\n")

    print("Original Uniform Grid Points:")
    # Using np.round to make the output cleaner and easier to read
    print(np.round(uniform_grid, 4))
    print("\nResulting Non-Uniform (Stretched) Grid Points:")
    print(np.round(stretched_grid, 4))
    print("\nNote how the points in the stretched grid are clustered near -1.0 and 1.0.")


# --- Parameters for the demonstration ---
# Number of points in our 1D grid
N = 15
# How much to cluster the points at the ends. Higher value = more clustering.
delta = 2.5

# Generate and print the grid points
generate_stretched_grid(N, delta)