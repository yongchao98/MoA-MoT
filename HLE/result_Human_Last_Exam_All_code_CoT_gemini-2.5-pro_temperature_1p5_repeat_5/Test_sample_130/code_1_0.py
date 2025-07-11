import numpy as np
from scipy.ndimage import median_filter

def print_image_with_header(header, img):
    """Helper function to print the image with a formatted header."""
    print(header)
    # The final code needs to output each number, so we print the array.
    for row in img:
        print(' '.join(map(str, row)))
    print("\n")

# --- Simulation Parameters ---
# Image size
N = 20
# Radius of the median filter's neighborhood
delta = 2
# The size of the square neighborhood is (2*delta + 1)
filter_size = 2 * delta + 1
# Total number of iterations to run
max_iterations = 5

# --- Initial Image Setup (t=0) ---
# Create a 20x20 black image (all zeros)
image = np.zeros((N, N), dtype=int)
# Add a 10x10 white square (all ones) in the center
image[5:15, 5:15] = 1

print_image_with_header("--- Image at t=0 ---", image)

# --- Iterative Filtering ---
current_image = image.copy()
for t in range(1, max_iterations + 1):
    # Apply the median filter.
    # The `size` parameter creates a square neighborhood of size (filter_size x filter_size).
    # This filter replaces each pixel with the median of its neighbors.
    current_image = median_filter(current_image, size=filter_size)

    # We print the image at specific time steps to observe the evolution.
    if t == 2 or t == max_iterations:
        header = f"--- Image at t={t} ---"
        print_image_with_header(header, current_image)
