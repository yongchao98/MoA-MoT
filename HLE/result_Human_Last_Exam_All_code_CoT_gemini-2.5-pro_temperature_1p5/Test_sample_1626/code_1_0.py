import numpy as np

def demonstrate_grid_spacing_with_tanh():
    """
    This function demonstrates how the tanh function can be used to
    control grid spacing by transforming a uniform grid into a
    non-uniform one with points clustered at the ends.
    """

    # 1. Create a uniform grid of 11 points from -3 to 3.
    # This represents our initial, uniformly-spaced logical coordinates.
    uniform_grid = np.linspace(-3, 3, 11)

    # 2. Apply the tanh function to the uniform grid.
    # This function smoothly maps the infinite real line to the interval (-1, 1).
    # This mapping is what 'stretches' the grid.
    transformed_grid = np.tanh(uniform_grid)

    # 3. Print the results to show the effect.
    print("This script demonstrates how the tanh function controls grid spacing.")
    print("-" * 60)
    print("Original Uniform Grid Points:")
    print(np.round(uniform_grid, 4))
    print("\nTransformed Grid Points using tanh:")
    print("Notice how the points are now clustered towards -1.0 and 1.0.")
    print(np.round(transformed_grid, 4))
    print("-" * 60)
    
    # Calculate and print the spacing between points to make it clearer
    uniform_spacing = np.diff(uniform_grid)
    transformed_spacing = np.diff(transformed_grid)
    
    print("\nSpacing between points in the Original Uniform Grid:")
    print(np.round(uniform_spacing, 4))
    print("\nSpacing between points in the Transformed Grid:")
    print("The spacing is small at the ends and large in the middle.")
    print(np.round(transformed_spacing, 4))


demonstrate_grid_spacing_with_tanh()