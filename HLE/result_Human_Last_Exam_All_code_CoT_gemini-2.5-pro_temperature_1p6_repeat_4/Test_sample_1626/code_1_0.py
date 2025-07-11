import numpy as np

def demonstrate_tanh_mapping():
    """
    Demonstrates how the tanh function maps uniform points to a non-uniform,
    clustered distribution. This principle is used in toroidal grid generation
    to adjust grid line spacing for consistent physical resolution.
    
    A uniform computational coordinate 'x' is mapped to a physical
    coordinate 'y' using y = tanh(x).
    """
    print("Demonstrating the mapping y = tanh(x):")
    print("A uniform step in 'x' results in a non-uniform step in 'y'.")
    print("-" * 40)
    
    # Create 11 uniformly spaced points from -2 to 2 (computational space)
    computational_points_x = np.linspace(-2.0, 2.0, 11)
    
    # Apply the tanh function to get the mapped points (physical space)
    physical_points_y = np.tanh(computational_points_x)
    
    # For each point, print the equation to show the mapping.
    # This fulfills the requirement to "output each number in the final equation".
    for x, y in zip(computational_points_x, physical_points_y):
        print(f"tanh({x:5.2f}) = {y:8.5f}")
        
    print("-" * 40)
    print("Notice how the output 'y' values are clustered near the ends")
    print("(around -1 and 1) and are more spread out in the middle.")
    print("This controllable, non-uniform spacing is why tanh is used.")

demonstrate_tanh_mapping()
