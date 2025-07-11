import numpy as np
from scipy.ndimage import median_filter

def print_image(img, title=""):
    """Helper function to print the image with a title."""
    print(title)
    # Print a text representation of the image array
    for row in img:
        print(" ".join(['■' if p == 1 else '□' for p in row]))
    print("\n" + "="*45 + "\n")

# --- 1. Setup the Simulation ---
# Image Size (N x N)
N = 25
# Neighborhood radius for the filter. We'll use a square of size (2*delta+1)
# which approximates a circle with radius delta.
delta = 2
filter_size = 2 * delta + 1

# --- 2. Create the Initial Image ---
# We create a binary image with a non-convex shape (a plus sign '+').
# This shape has sharp concave and convex corners, making it a good test case.
image = np.zeros((N, N), dtype=int)
center = N // 2
arm_length = N // 4
# Horizontal bar of the plus sign
image[center, center - arm_length : center + arm_length + 1] = 1
# Vertical bar of the plus sign
image[center - arm_length : center + arm_length + 1, center] = 1


# --- 3. Run the Simulation ---
# The core of the simulation is repeatedly applying the median filter.
# As we iterate (as t increases), the edges will change.

# Print the initial state (t=0)
print_image(image, "Initial Image (t = 0):")

# Apply the filter for a few iterations
current_image = image.copy()
num_iterations_to_show = [1, 3, 7]

for t in range(1, max(num_iterations_to_show) + 1):
    # Apply the local median filter. mode='constant' handles pixels near the border.
    current_image = median_filter(current_image, size=filter_size, mode='constant', cval=0)

    # Print the state at specific iterations
    if t in num_iterations_to_show:
        title = f"Image after t = {t} iterations:"
        print_image(current_image, title)

# --- 4. Final Explanation ---
# The code above simulates the process. As you can see from the outputs,
# the initial sharp '+' shape has its inner (concave) corners filled in and
# its outer (convex) corners rounded off. The shape becomes more compact
# and blob-like, tending towards a circle. This demonstrates how the edges
# are smoothed out as t increases.
