import numpy as np

def demonstrate_tanh_grid_spacing():
    """
    This function demonstrates how the hyperbolic tangent (tanh) function
    can be used to transform a uniformly spaced grid into a grid with
    points clustered at the boundaries. This is a common technique used
    to control grid resolution.
    """
    
    # 1. Create a set of 11 uniformly spaced points in a given domain (e.g., -3 to 3).
    # This represents the initial, logical grid space.
    uniform_points = np.linspace(-3, 3, 11)

    # 2. Apply the tanh function to these points.
    # The output will be in the range (-1, 1).
    tanh_transformed_points = np.tanh(uniform_points)
    
    # There is no single equation to solve, but we can display the transformation.
    # We will print each number involved in this demonstration.
    
    print("This script demonstrates the effect of the tanh function on grid spacing.")
    print("-" * 60)
    
    print("Original Uniformly Spaced Points:")
    # Printing each number from the original set
    for point in uniform_points:
        print(f"{point:9.4f}")
        
    print("\nTransformed Points using tanh(x):")
    print("Notice how the points are now clustered towards the ends (-1.0 and 1.0).")
    # Printing each number from the transformed set
    for point in tanh_transformed_points:
        print(f"{point:9.4f}")
    print("-" * 60)


if __name__ == "__main__":
    demonstrate_tanh_grid_spacing()
