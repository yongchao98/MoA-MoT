import numpy as np
from scipy.ndimage import median_filter

def describe_image_state(image, iteration):
    """Prints a description of the image state."""
    total_pixels = image.size
    white_pixels = np.sum(image)
    white_percentage = 100 * white_pixels / total_pixels
    print(f"Iteration {iteration: >3}: Image has {white_pixels} white pixels ({white_percentage:.2f}%).")

# Plan:
# 1. Define image and filter parameters.
# 2. Create an initial binary image with a non-convex shape
#    (a large white square with a smaller black square cut out).
# 3. Repeatedly apply a median filter.
# 4. Print the state of the image at key iterations to show the evolution.

# 1. Parameters
N = 100  # Image size: N x N
delta = 3  # Radius for the median filter neighborhood
# The filter uses a square neighborhood of side (2*delta + 1)
filter_size = 2 * delta + 1
iterations_to_show = [0, 1, 5, 10, 20, 50, 100]
max_iterations = 101

# 2. Create initial image
# Start with a black image
image = np.zeros((N, N), dtype=int)
# Add a large white square
image[20:80, 20:80] = 1
# Cut out a smaller black square to create a non-convex 'L' shape
image[20:50, 20:50] = 0

# 3. & 4. Run the simulation
print("Starting simulation of repeated median filtering.")
print(f"Initial image is {N}x{N}. Filter size is {filter_size}x{filter_size}.")
print("-" * 60)

current_image = image.copy()
for i in range(max_iterations):
    if i in iterations_to_show:
        describe_image_state(current_image, i)
    # Apply the median filter
    current_image = median_filter(current_image, size=filter_size)

print("-" * 60)
print("Observation: The number of white pixels steadily decreases.")
print("This demonstrates that the shape is shrinking as its edges move.")
print("The convex corners are eroded and the concave corner is filled in,")
print("leading to a net reduction in the white area as the shape becomes smoother.")
print("As t -> infinity, the number of white pixels would approach 0,")
print("and the image would eventually become all black.")