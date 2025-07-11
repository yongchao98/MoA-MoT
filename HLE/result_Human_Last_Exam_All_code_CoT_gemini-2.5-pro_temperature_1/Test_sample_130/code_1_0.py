import numpy as np
# This script requires the scipy library for the median_filter function.
# You can install it with: pip install scipy
from scipy.ndimage import median_filter

def print_image(img, title):
    """Prints a text representation of the binary image."""
    print(title)
    for row in img:
        # Use '#' for white pixels (1) and '.' for black pixels (0)
        print("".join(['#' if p == 1 else '.' for p in row]))
    print("\n")

# --- Simulation Parameters ---
# The local median function is defined as f(x) = Median over the neighborhood ||x-y|| <= delta.
# We will use the following parameter for delta:
delta = 2

# Image dimensions for the simulation
height, width = 20, 40

# We will show the image at these specific iteration numbers
iterations_to_show = [0, 5, 10, 15]

# --- Create Initial Image ---
# An initial image 'I' with a white cross shape on a black background.
image = np.zeros((height, width), dtype=int)
# Vertical bar of the cross
image[height//2 - 5 : height//2 + 5, width//2 - 1 : width//2 + 1] = 1
# Horizontal bar of the cross
image[height//2 - 1 : height//2 + 1, width//2 - 7 : width//2 + 7] = 1

# --- Define Median Filter Footprint ---
# A circular footprint is created to represent the neighborhood ||x-y|| <= delta.
y, x = np.ogrid[-delta:delta+1, -delta:delta+1]
footprint = x**2 + y**2 <= delta**2

# --- Run Simulation and Print Results ---
current_image = image.copy()
iteration_counter = 0

print("This script demonstrates the effect of an iterative local median filter on a binary image.")
print(f"The simulation uses a neighborhood radius (delta) of: {delta}\n")

while iteration_counter <= max(iterations_to_show):
    # Print the image at specified iterations
    if iteration_counter in iterations_to_show:
        print_image(current_image, f"Image at iteration t = {iteration_counter}:")

    # Stop the simulation if the shape has completely vanished
    if np.sum(current_image) == 0:
        print(f"The shape vanished completely after iteration {iteration_counter}.")
        break

    # Apply the median filter to get the image for the next iteration
    current_image = median_filter(current_image, footprint=footprint, mode='constant', cval=0)
    iteration_counter += 1
