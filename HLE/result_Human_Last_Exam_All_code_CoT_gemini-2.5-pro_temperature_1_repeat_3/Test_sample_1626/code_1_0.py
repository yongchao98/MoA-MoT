import numpy as np
import math

def demonstrate_grid_stretching():
    """
    This function demonstrates how the tanh function is used to control grid spacing.
    It creates a uniform grid and then transforms it into a non-uniform grid
    where points are clustered at the ends.
    """
    print("In toroidal grid generation, a stretching function is often used to map a uniform computational grid to a non-uniform physical grid.")
    print("This allows for higher resolution (denser grid points) in specific areas of interest, like near boundaries.")
    print("\nThe hyperbolic tangent function, tanh, is commonly used for this purpose.\n")

    # Parameters for the demonstration
    num_points = 15
    # A higher factor leads to more intense clustering at the ends
    stretching_factor = 3.0

    # 1. Create a uniform grid in the computational domain [-1, 1]
    uniform_grid = np.linspace(-1, 1, num_points)

    # 2. Apply a tanh-based transformation to create the non-uniform grid.
    # The formula `X_physical = tanh(k * X_comp) / tanh(k)` maps [-1, 1] to [-1, 1]
    # while clustering points towards the ends.
    denominator = math.tanh(stretching_factor)
    non_uniform_grid = np.tanh(stretching_factor * uniform_grid) / denominator

    print("--- Demonstration ---")
    print(f"Original uniform grid points between -1 and 1:")
    # Per the instructions, printing the numbers involved
    print([f"{x:.4f}" for x in uniform_grid])

    print(f"\nTransformed non-uniform grid points using tanh (stretching_factor={stretching_factor}):")
    # Per the instructions, printing the numbers involved
    print([f"{x:.4f}" for x in non_uniform_grid])
    print("\nNote how the points in the transformed grid are clustered near -1.0 and 1.0.")
    print("---------------------\n")


    # Final Answer
    # Identify the function from the provided choices.
    answer = "tanh"
    print("From the given choices, the specific mathematical function that governs grid spacing to ensure resolution consistency is:")
    print(answer)


if __name__ == '__main__':
    demonstrate_grid_stretching()

<<<C>>>