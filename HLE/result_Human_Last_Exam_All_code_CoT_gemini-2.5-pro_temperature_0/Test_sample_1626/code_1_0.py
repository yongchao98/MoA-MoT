import numpy as np

# The user is asking for the name of a mathematical function used for controlling
# grid spacing in toroidal grid generation. The correct answer is the hyperbolic
# tangent function, or 'tanh'.
#
# This script will demonstrate how the 'tanh' function can be used to "stretch"
# a grid, which is its primary role in this context. It takes a uniform set of
# points and redistributes them to be clustered towards the boundaries.

def demonstrate_grid_stretching(function_name, num_points=11):
    """
    Demonstrates how a function can create non-uniform grid spacing.

    Args:
        function_name (str): The name of the function to use.
        num_points (int): The number of grid points to generate.
    """
    print(f"Demonstrating grid spacing control using the '{function_name}' function.")

    # 1. Create a uniform grid in an abstract "computational space".
    # We use a range from -3 to 3 to clearly show the tanh compression effect.
    computational_grid = np.linspace(-3, 3, num_points)

    # 2. Apply the hyperbolic tangent function to the computational grid.
    # The result is a "physical grid" where points are clustered towards -1 and 1.
    if function_name == 'tanh':
        physical_grid = np.tanh(computational_grid)
    else:
        # Placeholder for other functions, not used here.
        physical_grid = computational_grid

    # 3. Print the results to show the effect.
    print("\nOriginal uniform grid (Computational Space):")
    # The equation for this grid is just a linear sequence.
    # We will print each number in the resulting array.
    for number in computational_grid:
        print(f"{number:6.2f}", end=" ")
    print()

    print(f"\nStretched grid using {function_name} (Physical Space):")
    # The equation used here is: y = tanh(x)
    # We will print each number from the final result.
    for number in physical_grid:
        print(f"{number:6.2f}", end=" ")
    print("\n\nAs you can see, the points in the second grid are clustered near the ends (-1.00 and 1.00).")
    print(f"This demonstrates how '{function_name}' governs grid spacing for resolution consistency.")


# The name of the function in question is 'tanh'.
correct_function_name = "tanh"

# Run the demonstration.
demonstrate_grid_stretching(correct_function_name)
