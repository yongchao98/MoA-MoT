import numpy as np
from scipy.ndimage import median_filter

def print_image(image, title):
    """Helper function to print the image array with a title."""
    print(title)
    # Replace 0s with '.' for better visual clarity
    for row in image:
        print(" ".join(['1' if pixel == 1 else '.' for pixel in row]))
    print("\n" + "="*40 + "\n")

# --- Configuration ---
# Image size
N = 20
# Size of the median filter neighborhood (e.g., 5x5). Must be odd.
# This corresponds to delta in the problem description.
filter_size = 5

# --- Initialization ---
# Create a black image (all zeros)
image = np.zeros((N, N), dtype=int)

# Create a 10x10 white square in the center
# The square has sharp 90-degree corners.
start, end = (N - 10) // 2, (N + 10) // 2
image[start:end, start:end] = 1

# --- Simulation ---
# 1. Print the initial state
print_image(image, "Initial Image (t=0): A 10x10 square.")

# 2. Apply the filter for a few iterations
num_iterations_1 = 4
for _ in range(num_iterations_1):
    image = median_filter(image, size=filter_size, mode='constant', cval=0)

print_image(image, f"After {num_iterations_1} iterations: Corners are rounded.")

# 3. Apply the filter for several more iterations
num_iterations_2 = 8
for _ in range(num_iterations_2):
    image = median_filter(image, size=filter_size, mode='constant', cval=0)

print_image(image, f"After {num_iterations_1 + num_iterations_2} total iterations: The shape has shrunk significantly.")
