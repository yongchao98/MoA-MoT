import numpy as np

def demonstrate_tanh_grid_spacing():
    """
    Demonstrates how the tanh function can be used to control grid spacing.
    A uniform grid in "computational space" is transformed into a
    non-uniform grid in "physical space" with points clustered at the center.
    """
    # 1. Create a uniform grid (e.g., 21 points from -4 to 4)
    # This represents our simple "computational" grid coordinates.
    computational_grid = np.linspace(-4, 4, 21)

    # 2. The core transformation using the tanh function.
    # This is the function that governs the new spacing.
    # y = tanh(x)
    physical_grid = np.tanh(computational_grid)

    print("This script demonstrates using the tanh function for grid spacing control.")
    print("-" * 60)
    print("Original Uniform Grid Points (Computational Space):")
    # The user instruction asked to "output each number in the final equation!"
    # We will show the input and output numbers of the transformation y = tanh(x).
    for point in computational_grid:
        print(f"{point:6.2f}", end=" ")
    print("\n")


    print("Transformed Non-Uniform Grid Points (Physical Space using tanh):")
    for point in physical_grid:
        print(f"{point:6.2f}", end=" ")
    print("\n")
    print("-" * 60)
    print("Notice how the points in the transformed grid are now clustered around 0.00,")
    print("achieving higher resolution in the center of the domain.")


if __name__ == '__main__':
    demonstrate_tanh_grid_spacing()