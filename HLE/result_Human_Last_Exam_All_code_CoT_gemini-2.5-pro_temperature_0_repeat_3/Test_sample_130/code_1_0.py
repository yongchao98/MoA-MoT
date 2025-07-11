import numpy as np
from scipy.ndimage import median_filter

def simulate_median_filter_evolution():
    """
    Simulates the long-term application of a local median filter on a binary image
    and explains the effect on the edges.
    """
    # --- Step 1: Define Parameters ---
    # These are the "numbers in the equation" f_t+delta(x)
    N = 100  # Image size: N x N
    delta = 3 # Radius of the median filter neighborhood, ||x-y|| <= delta

    # --- Step 2: Create an Initial Binary Image ---
    # We create a shape with both convex and concave corners (a plus sign)
    # to observe the evolution of its edges.
    image = np.zeros((N, N), dtype=np.uint8)
    center = N // 2
    width = N // 10
    length = N // 3
    # Create the vertical bar of the plus sign
    image[center - length : center + length, center - width : center + width] = 1
    # Create the horizontal bar of the plus sign
    image[center - width : center + width, center - length : center + length] = 1

    # --- Step 3: Define the Filter Footprint ---
    # The condition ||x-y|| <= delta defines a circular neighborhood.
    # We create a disk-shaped footprint for the median filter.
    y, x = np.ogrid[-delta:delta + 1, -delta:delta + 1]
    disk_footprint = x**2 + y**2 <= delta**2

    # --- Step 4: Explain and Run the Simulation ---
    print("--- Analyzing the operator f(I) = MedianFilter(I) over time t ---")
    print(f"Image Size (N): {N}")
    print(f"Filter Radius (delta): {delta}")
    print("The operator acts as a majority filter on the binary image.")
    print("We will track the area (number of white pixels) of the shape over iterations (t).")
    print("-" * 60)

    current_image = image.copy()
    initial_area = np.sum(current_image)
    print(f"t = 0 (Initial State): Area = {initial_area}")
    print("The initial shape is a plus sign, which has convex outer corners and concave inner corners.")

    # We run the simulation and print the area at specific timesteps t.
    max_iterations = 100
    iterations_to_show = [1, 5, 10, 25, 50, 75]
    
    for t in range(1, max_iterations + 1):
        # Apply the local median operator
        current_image = median_filter(current_image, footprint=disk_footprint)
        
        if t in iterations_to_show:
            area = np.sum(current_image)
            print(f"t = {t}: Area = {area}")
            if t == 1:
                print("  - Notice the sharp corners are immediately rounded.")
            if t == 10:
                print("  - The shape is becoming more circular as convex parts shrink and concave parts fill.")
            if area == 0:
                print("  - The shape has completely disappeared.")
                break
    
    print("-" * 60)
    print("--- Final Conclusion ---")
    print("The repeated application of the local median filter causes edges to move in a way that minimizes their length.")
    print("This process is known as 'motion by mean curvature'.")
    print("As t -> infinity, the final state of the edges depends on their initial topology:")
    print("1. Closed Edges (forming a shape): The boundary smooths out, and the shape shrinks until it disappears.")
    print("2. Open Edges (crossing the image): The boundary smooths out until it becomes a straight line.")

# Execute the simulation
simulate_median_filter_evolution()