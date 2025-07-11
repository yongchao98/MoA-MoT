import numpy as np

def demonstrate_tanh_grid_stretching():
    """
    This function demonstrates how the hyperbolic tangent (tanh) function
    is used to govern grid spacing. It transforms a uniform set of points
    into a non-uniform set, clustering them around the center. This is a
    key technique for ensuring resolution consistency in complex geometries
    like a torus.
    """
    print("--- Demonstrating Grid Spacing with tanh ---")

    # Define a set of 11 evenly spaced points in a computational domain.
    # In a real grid, this might represent an angular or radial index.
    uniform_points = np.linspace(-2.5, 2.5, 11)

    print("\n1. A set of uniform, linearly spaced points:")
    print(np.round(uniform_points, 4))
    print("\nThis uniform spacing can lead to resolution issues in a toroidal geometry.")

    # Apply the tanh function to transform the points.
    # This simulates the mapping to the physical grid coordinates.
    transformed_points = np.tanh(uniform_points)

    print("\n2. The same points transformed by the tanh function:")
    print(np.round(transformed_points, 4))
    print("\nNotice the points are now clustered around 0.0, with spacing increasing towards the edges.")

    # Show the "equation" for each point transformation.
    print("\n--- The tanh function governing each point's final position ---")
    for i in range(len(uniform_points)):
        original = uniform_points[i]
        transformed = transformed_points[i]
        # We print each number in the final equation: new_val = tanh(original_val)
        print(f"Point {i+1}: {transformed:8.4f} = tanh({original:5.2f})")

demonstrate_tanh_grid_stretching()