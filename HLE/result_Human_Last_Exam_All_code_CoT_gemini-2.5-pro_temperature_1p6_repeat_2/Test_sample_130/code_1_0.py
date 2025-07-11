import numpy as np
from scipy.ndimage import median_filter

def print_image(img, title):
    """Helper function to print the image array with a title and readable formatting."""
    print(title)
    # The final equation is the state of the image, represented by its pixel values.
    # We print each number in this final state.
    for row in img:
        print(" ".join(map(str, row)))
    print("\n" + "="*50 + "\n")

# --- Conceptual Explanation ---
# As we apply the local median filter iteratively, the edges of the image evolve in a process
# that smooths them out. This is known as motion by mean curvature.
# 1. Convex corners (like the tips of a plus sign) are eroded because a pixel at the tip
#    is surrounded by more 'background' pixels. Its median becomes 0, and it is removed.
# 2. Concave corners (like the inner corners of a plus sign) are filled in. A background
#    pixel in a concave nook is surrounded by more 'shape' pixels. Its median becomes 1,
#    and it becomes part of the shape.
# 3. Straight edges are stable and do not move, as a pixel on a straight edge has a balanced
#    neighborhood.
#
# The long-term effect (as t -> infinity) is that boundaries attempt to straighten out.
# Our 'plus' sign will first evolve into a square, and this square (a convex shape)
# will then shrink until it disappears, leaving a completely black image.

# --- Step 1: Create the initial image ---
N = 22
image = np.zeros((N, N), dtype=int)

# Create a "plus" sign shape to demonstrate the evolution of both
# convex and concave corners.
image[6:16, 9:13] = 1
image[9:13, 6:16] = 1

# --- Step 2 & 3: Apply the median filter iteratively ---
# The filter radius delta is 1, so the neighborhood size is (2*1+1)x(2*1+1) = 3x3.
filter_size = 3

# We will show the image at several time steps to visualize the evolution.
iterations_to_show = [0, 5, 15, 25]
current_image = image.copy()

print("This script demonstrates what happens when a local median filter is applied iteratively to a binary image.")
print("We start with a 'plus' shape and observe how its edges change over time.\n")

print_image(current_image, "Image at t=0 (Initial 'Plus' Shape):")

for t in range(1, max(iterations_to_show) + 1):
    # Apply the median filter
    current_image = median_filter(current_image, size=filter_size, mode='constant', cval=0)
    if t in iterations_to_show:
        # The equation for the image at this step is its matrix of pixel values.
        print_image(current_image, f"Final Equation of the Image at t={t}:")
