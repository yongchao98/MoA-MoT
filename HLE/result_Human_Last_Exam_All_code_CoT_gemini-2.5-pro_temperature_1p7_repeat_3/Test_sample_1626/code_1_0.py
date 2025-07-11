import math
import numpy as np

def demonstrate_tanh_spacing(num_points=21, packing_factor=2.0):
    """
    Demonstrates how the tanh function clusters grid points.

    Args:
        num_points (int): The number of grid points to generate.
        packing_factor (float): Controls the degree of clustering. Higher is more packed.
    """
    print(f"Demonstrating grid spacing with tanh function for {num_points} points.")
    print(f"Packing factor (alpha) = {packing_factor}\n")
    
    # 1. Create a uniform grid in the "computational space" from -1 to 1
    uniform_grid = np.linspace(-1, 1, num_points)
    
    # 2. Create the transformed grid in the "physical space"
    # The division by tanh(packing_factor) normalizes the output to the [-1, 1] range.
    tanh_normalizer = math.tanh(packing_factor)
    transformed_grid = np.tanh(packing_factor * uniform_grid) / tanh_normalizer
    
    # 3. Print the results for comparison
    print(f"{'Uniform Grid':<25} | {'Tanh-Spaced Grid':<25}")
    print("-" * 55)
    for i in range(num_points):
        print(f"{uniform_grid[i]:<25.8f} | {transformed_grid[i]:<25.8f}")
        
    print("\nNotice how the values in the 'Tanh-Spaced Grid' are clustered near -1.0 and 1.0,")
    print("providing higher resolution at the boundaries of the domain.")

if __name__ == '__main__':
    demonstrate_tanh_spacing()
