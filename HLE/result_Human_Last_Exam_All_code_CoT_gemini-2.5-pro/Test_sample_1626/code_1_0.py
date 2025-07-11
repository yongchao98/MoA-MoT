import numpy as np

def demonstrate_tanh_stretching():
    """
    Demonstrates how the tanh function is used to create a non-uniform grid
    by clustering points around the center.
    """
    # Number of grid points
    num_points = 21

    # A stretching factor. Higher values lead to more clustering.
    alpha = 2.5

    # 1. Create a uniform grid in "logical space" from -1 to 1
    # This represents our initial, evenly spaced coordinates.
    uniform_grid = np.linspace(-1.0, 1.0, num_points)

    # 2. Apply the tanh transformation to create the "physical space" grid
    # The division by np.tanh(alpha) normalizes the output to remain in the [-1, 1] range.
    stretched_grid = np.tanh(alpha * uniform_grid) / np.tanh(alpha)

    print("--- Demonstration of Grid Stretching with tanh ---")
    print(f"Stretching factor (alpha): {alpha}\n")

    print("Original Uniform Grid vs. Tanh-Stretched Grid:")
    print("-" * 50)
    print(f"{'Index':<6} | {'Uniform Point':<15} | {'Stretched Point':<20}")
    print("-" * 50)

    for i in range(num_points):
        # The prompt asks to output each number in the final equation.
        # Here we print the transformation for each point.
        # Equation: stretched_grid[i] = tanh(alpha * uniform_grid[i]) / tanh(alpha)
        
        u_val = uniform_grid[i]
        s_val = stretched_grid[i]
        
        # Printing the components of the calculation for each point
        print(f"i={i:<4} | {u_val:<15.4f} | {s_val:<20.4f}")
        
    print("-" * 50)
    print("\nNotice how the 'Stretched Point' values are clustered much more")
    print("densely around 0.0 compared to the 'Uniform Point' values.")


if __name__ == "__main__":
    demonstrate_tanh_stretching()
