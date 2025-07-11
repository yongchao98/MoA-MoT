import numpy as np

def demonstrate_tanh_stretching(num_points=11, stretching_factor=2.0):
    """
    Demonstrates how the tanh function is used to create a non-uniform grid
    by stretching a uniform grid.

    Args:
        num_points (int): The number of grid points to generate.
        stretching_factor (float): A factor controlling the degree of clustering.
                                   Higher values lead to more clustering at the ends.
    """
    # 1. Create a uniform grid in "computational space" from -1 to 1.
    # This is our starting point with perfectly consistent spacing.
    uniform_grid = np.linspace(-1.0, 1.0, num_points)

    # 2. Apply the tanh mapping function to get the "physical space" grid.
    # The equation used is a common formulation for grid stretching.
    # The denominator normalizes the output to ensure the new grid also spans [-1, 1].
    numerator = np.tanh(stretching_factor * uniform_grid)
    denominator = np.tanh(stretching_factor)
    stretched_grid = numerator / denominator

    # The equation is: stretched_coord = tanh(stretching_factor * uniform_coord) / tanh(stretching_factor)
    # Let's print the numbers that go into the equation for the first and last points as an example.
    
    print("--- Grid Spacing Demonstration using tanh ---")
    print(f"\nEquation form: new_x = tanh(factor * old_x) / tanh(factor)")
    print(f"Stretching Factor: {stretching_factor}\n")

    print(f"{'Original Uniform Point':<25} | {'Transformed Stretched Point':<25}")
    print("-" * 55)
    for i in range(num_points):
        # We demonstrate the calculation for each point as requested
        # For each point, the 'equation' is:
        # stretched_grid[i] = tanh(stretching_factor * uniform_grid[i]) / tanh(stretching_factor)
        # We will output the numbers involved in this equation:
        # uniform_grid[i], stretching_factor, and stretched_grid[i]
        
        print(f"Input: {uniform_grid[i]:<19.4f} | Output: {stretched_grid[i]:<25.4f}")

    print("\nNote how the output points are clustered near -1.0 and 1.0,")
    print("while the points near the center are more spread out.")


if __name__ == '__main__':
    demonstrate_tanh_stretching()