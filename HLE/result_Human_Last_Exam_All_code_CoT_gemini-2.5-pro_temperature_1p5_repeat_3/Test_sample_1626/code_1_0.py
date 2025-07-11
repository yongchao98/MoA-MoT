import numpy as np

def generate_tanh_grid():
    """
    Demonstrates the use of the tanh function for grid spacing.
    This function generates a uniform grid and a tanh-stretched grid
    to show how points can be clustered in a region of interest.
    """
    # --- Parameters ---
    # Number of grid points
    num_points = 21
    # Half-width of the physical domain [-L, L]
    domain_half_width = 10.0
    # Clustering factor. Higher alpha = more clustering at the center.
    alpha = 2.5

    # 1. Create a uniform computational grid from -1 to 1
    computational_grid = np.linspace(-1, 1, num_points)

    # 2. Apply the tanh mapping function to get the physical grid
    # The equation is: x = L * tanh(alpha * xi) / tanh(alpha)
    # This maps the uniform 'computational_grid' (xi) to a non-uniform 'physical_grid' (x)
    tanh_denominator = np.tanh(alpha)
    physical_grid_tanh = domain_half_width * np.tanh(alpha * computational_grid) / tanh_denominator

    # 3. For comparison, create a perfectly uniform physical grid
    physical_grid_uniform = np.linspace(-domain_half_width, domain_half_width, num_points)

    # 4. Print the results
    # The "final equation" in this case is the transformation that produces the grid.
    # We will print the numbers that result from this equation: the grid points themselves.
    print(f"Demonstration of tanh grid spacing for a domain from {-domain_half_width} to {domain_half_width}")
    print("-" * 70)
    print("A. Uniformly spaced grid points:")
    print(np.round(physical_grid_uniform, 3))
    print("\nB. Tanh-stretched grid points (note the clustering around 0.0):")
    print(np.round(physical_grid_tanh, 3))
    print("-" * 70)

if __name__ == '__main__':
    generate_tanh_grid()
