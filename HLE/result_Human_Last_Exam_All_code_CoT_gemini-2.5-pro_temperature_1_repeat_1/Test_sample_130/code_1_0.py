import numpy as np
from scipy.ndimage import median_filter

def print_image(img, title):
    """
    Helper function to print the image with a title.
    Uses '█' for white pixels (1) and ' ' for black pixels (0) for better visual clarity.
    """
    print(title)
    for row in img:
        print("".join(['█' if pixel == 1 else ' ' for pixel in row]))
    print("\n" + "="*len(img[0]) + "\n")

# Plan:
# 1. Create a 20x20 binary image with an L-shape. This shape is ideal because it has
#    both convex (outer) and concave (inner) corners.
# 2. Define the neighborhood for the median filter. A size of 3x3 corresponds to delta=1.
# 3. Iteratively apply the median filter to the image.
# 4. Print the image at different stages (t=0, 1, 5, 10, 20) to show the evolution of the edges.

# Step 1: Create an initial binary image (I)
N = 20
image = np.zeros((N, N), dtype=int)
# Create the L-shape by setting pixel values to 1
image[5:15, 5:10] = 1  # Vertical part of the 'L'
image[10:15, 5:15] = 1 # Horizontal part of the 'L'

# Step 2: Define filter parameters
# A neighborhood size of 3 means a 3x3 window, which corresponds to delta=1.
neighborhood_size = 3

# --- Simulation ---

# Print the initial state of the image
print_image(image, "Image at t=0:")

# Step 3 & 4: Apply the filter iteratively and print the results
current_image = image.copy()
last_iteration_printed = 0

# We will show the image after 1, 5, 10, and 20 iterations.
for t in [1, 5, 10, 20]:
    # Apply the filter for the required number of new iterations
    num_steps = t - last_iteration_printed
    for _ in range(num_steps):
        # The median filter is applied here. mode='constant' handles pixels at the boundary
        # of the image by padding with a constant value (0, or black).
        current_image = median_filter(current_image, size=neighborhood_size, mode='constant', cval=0)
    
    # Print the resulting image
    print_image(current_image, f"Image after t={t} iterations:")
    last_iteration_printed = t
