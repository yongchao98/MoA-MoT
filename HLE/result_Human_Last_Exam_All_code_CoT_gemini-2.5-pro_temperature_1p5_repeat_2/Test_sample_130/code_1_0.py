import numpy as np
from scipy.ndimage import median_filter
import time
import os

def print_image(image, title):
    """Helper function to print the binary image to the console."""
    print(title)
    for row in image:
        print(" ".join(['■' if pixel == 1 else '·' for pixel in row]))
    print("\n" + "="*30 + "\n")

# --- Parameters ---
N = 20  # Image size
delta = 1  # Radius of the median filter neighborhood
filter_size = 2 * delta + 1
iterations_to_show = [0, 1, 5, 10, 15, 20]

# --- Initial Image Setup ---
# Create an empty image (all 0s)
image = np.zeros((N, N), dtype=int)

# Create an 'L' shape in the image
image[5:14, 5] = 1
image[13, 5:10] = 1

# --- Simulation ---
current_image = image.copy()
print_image(current_image, "Initial Image (t=0):")

for t in range(1, max(iterations_to_show) + 1):
    # Apply the local median filter
    current_image = median_filter(current_image, size=filter_size, mode='constant', cval=0)

    # Print the image at specified iterations
    if t in iterations_to_show:
        print_image(current_image, f"Image after {t} iterations:")
