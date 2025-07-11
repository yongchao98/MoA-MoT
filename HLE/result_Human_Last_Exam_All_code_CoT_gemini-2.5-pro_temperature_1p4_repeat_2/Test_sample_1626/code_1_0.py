import numpy as np

def demonstrate_tanh_grid_spacing():
    """
    Demonstrates how the tanh function is used to create a non-uniform grid
    by "packing" points around the center of the domain.
    """
    # 1. Define a uniform grid in "computational space" from -1 to 1.
    num_points = 11
    uniform_grid = np.linspace(-1.0, 1.0, num_points)

    # 2. Define the packing parameter. Higher values mean more packing.
    packing_factor = 2.5

    # 3. Apply the tanh transformation function.
    # The division by np.tanh(packing_factor) ensures the domain remains [-1, 1].
    # This is a common form of a grid stretching function.
    # The equation is: x_physical = tanh(packing_factor * x_computational) / tanh(packing_factor)
    denominator = np.tanh(packing_factor)
    stretched_grid = np.tanh(packing_factor * uniform_grid) / denominator

    print("This script demonstrates how a tanh function creates a non-uniform grid.")
    print("-" * 50)
    print(f"Using a packing factor of: {packing_factor}")
    print(f"The denominator tanh({packing_factor}) is: {denominator:.4f}")
    print("-" * 50)
    print("Original Uniform Grid Points:")
    for point in uniform_grid:
        print(f"{point:8.4f}")

    print("\nStretched Grid Points (governed by tanh):")
    # As requested, outputting each number in the final transformation.
    # Notice how the points are now clustered around 0.0.
    for point in stretched_grid:
        print(f"{point:8.4f}")

    print("\nObserve that the spacing between points in the 'Stretched Grid'")
    print("is smaller near the center (0.0) and larger near the ends (-1.0 and 1.0).")
    print("This is how tanh governs grid spacing to control resolution.")

if __name__ == '__main__':
    demonstrate_tanh_grid_spacing()