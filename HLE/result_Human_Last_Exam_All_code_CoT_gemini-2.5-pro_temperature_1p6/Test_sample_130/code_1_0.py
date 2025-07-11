import numpy as np
from scipy.ndimage import median_filter
import time

def print_image(image):
    """Prints the image with 1s as '#' and 0s as '.' for better visualization."""
    for row in image:
        print("".join(['#' if p == 1 else '.' for p in row]))
    print("\n")

def run_simulation():
    """
    Simulates the iterative application of a median filter on a binary image.
    """
    # We use a delta of 1, which means a 3x3 neighborhood for the median filter.
    # The operator is f_{t+δ}(x) = Median_{||x-y||<=δ}. With δ=1, this is a 3x3 median filter.
    delta = 1
    filter_size = 2 * delta + 1
    
    # Initialize a 16x16 image with a U-shape.
    # This shape has both convex (outer) and concave (inner) corners.
    image = np.zeros((16, 16), dtype=int)
    image[4:12, 4] = 1
    image[4:12, 11] = 1
    image[11, 4:12] = 1
    
    print("This script simulates the evolution of a binary image under an iterative median filter.")
    print(f"The filter uses a neighborhood radius delta = {delta}, resulting in a {filter_size}x{filter_size} filter.\n")
    
    print("Initial image (t=0):")
    # To satisfy the "output each number" requirement, we print the raw numpy array.
    # We also print a more visual version.
    print(image)
    print_image(image)
    
    # We iterate up to 20 times to see the evolution.
    max_iterations = 20
    for t in range(1, max_iterations + 1):
        # Store the previous image to check for convergence
        previous_image = np.copy(image)
        
        # Apply the median filter
        image = median_filter(image, size=filter_size, mode='constant', cval=0)
        
        print(f"Image at t={t}:")
        print(image)
        print_image(image)
        
        # Check if the image has stopped changing (reached a stable state)
        if np.array_equal(image, previous_image):
            print(f"The image has reached a stable state at t={t}.")
            break
            
    if t == max_iterations:
        print("Simulation finished after maximum iterations.")

if __name__ == '__main__':
    run_simulation()