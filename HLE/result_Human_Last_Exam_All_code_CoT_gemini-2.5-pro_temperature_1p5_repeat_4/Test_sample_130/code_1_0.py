import numpy as np
from scipy.ndimage import generic_filter

def simulate_median_filter_on_image():
    """
    Simulates the iterative application of a local median filter on a binary image.
    """
    # 1. Define image parameters
    N = 100  # Image size: N x N
    delta = 3  # Radius of the local median operator
    num_iterations = 50

    # 2. Create the initial binary image
    # A black image with a white square in the center
    image = np.zeros((N, N), dtype=np.uint8)
    square_size = 40
    start = (N - square_size) // 2
    end = start + square_size
    image[start:end, start:end] = 1

    # 3. Define the filter's neighborhood (footprint)
    # A circular neighborhood is defined by ||x-y|| <= delta
    y, x = np.ogrid[-delta:delta + 1, -delta:delta + 1]
    footprint = x**2 + y**2 <= delta**2

    # 4. Define the local median function for a binary image
    # For binary {0, 1} data, the median is 1 if there are more 1s than 0s,
    # and 0 otherwise. This is a majority filter.
    def majority_vote(neighborhood_1d):
        # Tie-breaking: if counts are equal, default to 0 (black).
        # This creates an "erosive" bias.
        ones = np.sum(neighborhood_1d)
        zeros = len(neighborhood_1d) - ones
        if ones > zeros:
            return 1
        else:
            return 0

    print("--- Simulating Iterative Local Median Filter ---")
    print(f"Initial Image: {N}x{N}, White Square Size: {square_size}x{square_size}")
    print(f"Filter Radius (delta): {delta}\n")
    print("The 'Equation' is I_t+1 = f(I_t), where f is the median filter.")
    print("We will print the total number of white pixels (Area) for I_t at each step t.\n")

    current_image = image.copy()

    # 5. Apply the filter iteratively
    for t in range(num_iterations + 1):
        # Calculate the area (sum of all pixel values)
        area = np.sum(current_image)
        
        # Output the result for the current step (the "number in the final equation")
        print(f"t = {t:2d}, Area (Number of white pixels) = {int(area)}")
        
        # Stop if the shape has disappeared
        if area == 0:
            print("\nThe white square has completely disappeared.")
            break
            
        # Apply the filter to get the next state of the image
        current_image = generic_filter(current_image, majority_vote, footprint=footprint)
        
    if np.sum(current_image) > 0:
        print("\nSimulation finished. Note how the area consistently decreases.")

simulate_median_filter_on_image()