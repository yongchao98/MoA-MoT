import numpy as np
from scipy.ndimage import median_filter

def print_image(img, title):
    """Helper function to print the image array in a readable format."""
    print(title)
    # Use '■' for white pixels and ' ' for black pixels for better visualization
    for row in img:
        print(" ".join(['■' if p == 1 else '·' for p in row]))
    print("-" * (2 * len(img[0]) - 1))

# 1. Define image and filter parameters
N = 21  # Image size (odd number to have a perfect center)
delta = 1  # Radius of the neighborhood, leading to a (2*delta+1)x(2*delta+1) filter
num_iterations = 10

# 2. Create an initial binary image
# A '+' shape is excellent for showing how convex and concave corners evolve.
image = np.zeros((N, N), dtype=int)
center = N // 2
arm_length = 5
# Horizontal arm
image[center, center - arm_length : center + arm_length + 1] = 1
# Vertical arm
image[center - arm_length : center + arm_length + 1, center] = 1

# 3. Print the initial state
print_image(image, "Initial Image (t=0):")

# 4. Apply the median filter iteratively
current_image = image.copy()
iterations_to_print = [1, 2, 5, 10]

for t in range(1, num_iterations + 1):
    # Apply the local median filter.
    # The 'size' parameter corresponds to the side length of the square neighborhood.
    # 'mode' defines how to handle boundaries; 'constant' fills with a constant value.
    current_image = median_filter(current_image, size=2*delta + 1, mode='constant', cval=0)

    # Print the image at specified iterations
    if t in iterations_to_print:
        # The output below shows the numbers in the final array for this step
        print_image(current_image, f"Image after {t} iterations:")
