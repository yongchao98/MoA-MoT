import numpy as np

# A function to demonstrate tanh stretching.
# It maps points from a uniform computational space 'xi' (from -c to c)
# to a physical space 'x' (from -1 to 1), clustering points around 0.
def tanh_stretch(xi, c):
    return np.tanh(c * xi) / np.tanh(c)

def main():
    # Number of grid points
    num_points = 11
    
    # A stretching factor. Higher values result in more clustering.
    stretching_factor = 2.0
    
    # 1. Create a uniform grid in the computational space [-1, 1]
    uniform_grid = np.linspace(-1.0, 1.0, num_points)
    
    print("--- Uniform Grid ---")
    print("This grid has evenly spaced points:")
    for point in uniform_grid:
        print(f"{point:8.4f}")
    
    print("\n--- Tanh-Stretched Grid ---")
    print("This grid, governed by the tanh function, clusters points near the center (0.0):")
    # 2. Apply the tanh stretching function to the uniform grid
    stretched_grid = tanh_stretch(uniform_grid, stretching_factor)
    for point in stretched_grid:
        print(f"{point:8.4f}")

if __name__ == "__main__":
    main()
