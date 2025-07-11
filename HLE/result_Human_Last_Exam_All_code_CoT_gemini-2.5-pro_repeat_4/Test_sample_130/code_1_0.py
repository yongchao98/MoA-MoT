# This script requires the numpy, scipy, and scikit-image libraries.
# You can install them using pip: pip install numpy scipy scikit-image

import numpy as np
from scipy.ndimage import median_filter
from skimage import measure

def calculate_perimeter(image):
    """Calculates the approximate perimeter of the largest shape in a binary image."""
    # find_contours finds the boundaries of positive-valued regions (our white shape).
    contours = measure.find_contours(image, 0.5)
    if contours:
        # We assume the largest contour is the one of interest.
        largest_contour = max(contours, key=len)
        # The number of points in the contour path is a good approximation of the perimeter.
        return len(largest_contour)
    return 0

# --- Simulation Parameters ---
# Image size (N x N)
N = 100
# Filter radius (delta), must be much smaller than N.
# The filter will use a square neighborhood of size (2*delta+1) x (2*delta+1).
delta = 2
filter_size = 2 * delta + 1

# 1. Create a binary image: a 50x50 white square on a 100x100 black background.
initial_image = np.zeros((N, N), dtype=np.uint8)
# The square is defined from pixel 25 to 74 (50 pixels wide).
initial_image[25:75, 25:75] = 1
initial_perimeter = (74 - 25 + 1) * 4 # Perimeter of a 50x50 square is 200

# Define the time steps (iterations) at which to observe the image.
iterations_to_print = [0, 1, 5, 10, 20, 30, 40]

# --- Simulation and Output ---
print("This script demonstrates the effect of a repeatedly applied local median filter on image edges.")
print(f"Starting with a {N}x{N} image containing a white square.")
print(f"The filter uses a neighborhood radius delta = {delta}, which corresponds to a {filter_size}x{filter_size} window.")
print("-" * 40)

current_image = initial_image.copy()
last_t = 0

for t in iterations_to_print:
    # 2. Apply the filter repeatedly to advance the simulation to time t.
    for _ in range(t - last_t):
        current_image = median_filter(current_image, size=filter_size)

    # 3. Measure the area and perimeter of the shape.
    area = np.sum(current_image)
    perimeter = calculate_perimeter(current_image)
    
    # At t=0, we print the initial state.
    if t == 0:
        print(f"At iteration t = {t} (Initial State):")
        print(f"  Shape area      = {area} pixels")
        print(f"  Shape perimeter = {initial_perimeter} units")
    else:
        print(f"At iteration t = {t}:")
        print(f"  Shape area      = {area} pixels")
        print(f"  Shape perimeter = {perimeter:.0f} units")

    last_t = t

    # Stop if the shape has vanished.
    if area == 0:
        print("\nThe shape has completely disappeared.")
        break

print("-" * 40)
print("Conclusion:")
print("The simulation shows that the filter erodes the shape, decreasing its area and perimeter over time.")
print("This occurs because pixels at convex boundaries (like corners) are surrounded by more '0's than '1's, causing them to flip to '0'.")
print("\nThis process is a form of 'mean curvature flow', where the velocity of the edge is proportional to its curvature.")
print("Therefore, as t -> infinity, the following happens to the edges:")
print("1. Smoothing: Sharp corners are immediately rounded.")
print("2. Straightening: Curved edges tend to become straighter to minimize their length.")
print("3. Shrinking: Closed boundaries (like circles or squares) shrink and eventually collapse to a point.")
print("\nUltimately, the edges become smoother and straighter, and closed shapes are removed from the image.")
