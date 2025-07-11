import numpy as np
from scipy.ndimage import median_filter
import sys

# This script demonstrates the effect of an iterative median filter on a binary image.
# We will observe the evolution of the edges of an L-shaped figure to understand
# the long-term behavior of this operation.

def print_image(img, title, stream=sys.stdout):
    """
    Prints a binary image to the console using block characters for better visualization.
    """
    stream.write(title + "\n")
    for row in img:
        # We use '██' for white pixels (1) and '  ' for black pixels (0).
        stream.write("".join(['██' if pixel == 1 else '  ' for pixel in row]) + "\n")
    stream.write("\n" + "="*41 + "\n\n")

# --- Parameters ---
# Image size
N = 20
# The operator radius, delta (δ). A delta of 1 means a 3x3 neighborhood.
delta = 1
# The filter size is (2*delta + 1) x (2*delta + 1)
filter_size = 2 * delta + 1
# We will show the image at these specific time steps (iterations)
iterations_to_show = [0, 1, 3, 10]

# --- Initial Image Setup (t=0) ---
# We create a binary image with a sharp L-shape. This allows us to see the
# effect on straight edges, a convex corner, and a concave corner simultaneously.
initial_image = np.zeros((N, N), dtype=int)
initial_image[5:15, 5:8] = 1   # Vertical bar of the 'L'
initial_image[12:15, 5:15] = 1 # Horizontal bar of the 'L'

current_image = initial_image.copy()

# --- Simulation Loop ---
# We now apply the median filter operator iteratively. With each step, the image
# I_t is transformed into I_{t+1}.
print(f"Applying a {filter_size}x{filter_size} median filter iteratively to an L-shaped figure.\n")

for t in range(max(iterations_to_show) + 1):
    # Check if the current iteration t is one we want to display
    if t in iterations_to_show:
        print_image(current_image, f"Image at iteration t = {t}")

    # This line applies the local median operator for one time step:
    # I_{t+1}(x) = Median_{||x-y||<=delta} (I_t(y))
    current_image = median_filter(current_image, size=filter_size, mode='constant', cval=0)

# The output above demonstrates the process.
# At t=1, the sharp corners begin to round.
# By t=3, the rounding is more pronounced, and the shape has started to shrink.
# By t=10, the shape is significantly smaller and smoother, having lost its initial sharp features.
# This illustrates the principle of curvature-driven flow, where edges move to become smoother.