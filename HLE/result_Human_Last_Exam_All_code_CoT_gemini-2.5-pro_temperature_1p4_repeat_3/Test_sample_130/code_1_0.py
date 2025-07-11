import numpy as np
from scipy.ndimage import median_filter

def print_image(image, title):
    """Helper function to print the image grid."""
    print(title)
    for row in image:
        # Represent 0 as '.' and 1 as '#' for better visibility
        print(" ".join(['#' if p == 1 else '.' for p in row]))
    print("-" * (2 * image.shape[1]))

# --- Parameters ---
# The size of the image grid
image_size = 21
# The radius of the median filter neighborhood. delta=1 means a 3x3 window.
delta = 1
# The number of times to apply the filter
iterations = 15

# --- Create the Initial Image ---
# An image with a cross shape, which has both convex and concave corners.
image = np.zeros((image_size, image_size), dtype=int)
center = image_size // 2
arm_length = 6
# Vertical arm
image[center - arm_length : center + arm_length + 1, center - 1 : center + 2] = 1
# Horizontal arm
image[center - 1 : center + 2, center - arm_length : center + arm_length + 1] = 1

# --- Print the Initial State ---
print_image(image, "Initial Image (t=0):")

# --- Apply the Filter Iteratively ---
processed_image = np.copy(image)
# The filter size is (2*delta + 1) x (2*delta + 1)
filter_size = 2 * delta + 1
for i in range(iterations):
    processed_image = median_filter(processed_image, size=filter_size, mode='constant', cval=0)
    # Optional: print intermediate steps to see the evolution
    # if (i + 1) % 5 == 0:
    #     print_image(processed_image, f"Image after {i+1} iterations:")


# --- Print the Final State ---
print_image(processed_image, f"Image after {iterations} iterations:")
