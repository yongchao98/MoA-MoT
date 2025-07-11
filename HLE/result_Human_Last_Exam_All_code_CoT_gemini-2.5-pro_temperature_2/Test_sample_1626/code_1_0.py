import math
import numpy as np

def demonstrate_tanh_for_grid_spacing():
    """
    This function demonstrates why the hyperbolic tangent (tanh) is an excellent
    choice for controlling grid spacing. It takes a uniformly spaced set of
    points and 'stretches' them, clustering the output points near the
    boundaries of the -1 to 1 interval.
    """

    print("--- Demonstrating Grid Point Clustering with tanh ---")
    print("In grid generation, we often start with uniform points in a 'computational space'")
    print("and map them to a 'physical space' where spacing is non-uniform.")
    print("\nThe function y = tanh(x) is used for this mapping.\n")
    
    # Generate 11 uniform points in our "computational space" from -2.5 to 2.5
    computational_points = np.linspace(-2.5, 2.5, 11)
    
    print("We will use the following equation to transform the points:")
    print("Physical_Point = tanh(Computational_Point)\n")

    print(f"{'Computational Point (x)':>25} | {'Physical Point (y = tanh(x))':>30}")
    print("-" * 60)
    
    # Calculate and print the mapping for each point
    # This fulfills the prompt to output each number in the final equation
    for x in computational_points:
        y = math.tanh(x)
        print(f"{x:>25.2f} | {y:>30.6f}")
        
    print("-" * 60)
    print("\nNotice how the output 'Physical Points' are clustered near -1.0 and 1.0,")
    print("while the points near the center are more spread out. This is the desired effect.")

demonstrate_tanh_for_grid_spacing()