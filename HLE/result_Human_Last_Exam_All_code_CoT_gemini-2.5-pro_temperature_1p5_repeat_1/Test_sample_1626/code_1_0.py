import numpy as np

def generate_stretched_grid(num_points, stretching_factor):
    """
    Generates a 1D grid stretched using the tanh function.

    Args:
        num_points (int): The number of grid points to generate.
        stretching_factor (float): A factor to control the degree of clustering.
                                   Higher values lead to more clustering.
    
    Returns:
        numpy.ndarray: An array of non-uniformly spaced grid points from -1.0 to 1.0.
    """
    # 1. Create a uniform grid in computational space (-1 to 1)
    uniform_grid = np.linspace(-1, 1, num_points)
    
    # 2. Apply the tanh stretching function
    # The denominator normalizes the output to ensure the domain remains [-1, 1]
    stretched_grid = np.tanh(stretching_factor * uniform_grid) / np.tanh(stretching_factor)
    
    return stretched_grid

# --- Parameters ---
# The total number of points in our grid
N_POINTS = 21 
# This parameter controls how much the grid is "bunched up" in the middle.
# A value of 0 would result in a uniform grid. Higher values increase clustering.
ALPHA = 2.0 

# --- Generate and Print the Grid ---
final_grid_points = generate_stretched_grid(N_POINTS, ALPHA)

print(f"Generated {N_POINTS} grid points using a tanh stretching function with a factor of {ALPHA}.")
print("The points are clustered around 0.0, providing higher resolution in the center:")
# The problem requests to output each number from the result
for point in final_grid_points:
    # We print each resulting grid point coordinate
    print(f"{point:.6f}")
