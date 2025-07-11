import numpy as np

def demonstrate_grid_spacing():
    """
    Demonstrates the use of the tanh function for grid spacing control.
    It generates and prints a linearly spaced grid and a tanh-spaced grid.
    """
    # Number of points in the grid
    num_points = 15
    
    # Stretching factor for the tanh function. A larger value means more clustering.
    stretching_factor = 2.5

    # 1. Generate a uniform computational grid from -1 to 1
    uniform_grid = np.linspace(-1, 1, num_points)

    # 2. Generate a stretched physical grid using the tanh function
    # The formula maps the uniform grid to a new grid where points are
    # clustered towards the ends (-1 and 1).
    numerator = np.tanh(stretching_factor * uniform_grid)
    denominator = np.tanh(stretching_factor)
    tanh_spaced_grid = numerator / denominator

    print("Demonstration of Grid Spacing Control:")
    print("-" * 40)
    
    print(f"A. Linearly Spaced Grid (Uniform Resolution):")
    # The equation for this is x_i = -1 + i * (2 / (N-1))
    # We will print each resulting number from this equation.
    print(np.round(uniform_grid, 4))
    print("\nNote how the distance between consecutive points is constant.")
    
    print("-" * 40)
    
    print(f"B. Tanh-Spaced Grid (Clustered Resolution):")
    # The equation is x_i = tanh(c * u_i) / tanh(c), where u is the uniform grid and c is the stretching factor.
    # We will print each resulting number from this equation.
    print(f"Stretching factor used in equation: {stretching_factor}")
    print(np.round(tanh_spaced_grid, 4))
    print("\nNote how the points are now clustered near the ends (-1.0 and 1.0).")
    print("-" * 40)

demonstrate_grid_spacing()