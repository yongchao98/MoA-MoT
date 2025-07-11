import numpy as np
from scipy.signal import convolve2d

def run_median_filter_simulation(image, delta, iterations):
    """
    Applies a local median (majority) filter to a binary image iteratively.

    Args:
        image (np.array): The initial binary image (0s and 1s).
        delta (int): The radius of the local neighborhood (e.g., delta=1 gives a 3x3 neighborhood).
        iterations (int): The number of times to apply the filter.

    Returns:
        np.array: The image after the final iteration.
    """
    # On a binary image, the median filter is a majority filter.
    # The neighborhood is a square of size (2*delta + 1) x (2*delta + 1).
    kernel_size = 2 * delta + 1
    kernel = np.ones((kernel_size, kernel_size))
    
    # The new pixel value is 1 if the sum of its neighbors (including itself)
    # is greater than half the number of pixels in the neighborhood.
    majority_threshold = kernel.size / 2.0
    
    img = image.copy()
    
    print("Initial Image:")
    # Print the image using '#' for 1s and '.' for 0s for clarity
    for row in img:
        print("".join(map(str, row)).replace('0', '.').replace('1', '#'))

    for i in range(iterations):
        # We use convolution to efficiently calculate the sum of neighbors for each pixel.
        # 'mode=same' ensures the output image has the same size as the input.
        # 'boundary=fill' handles edges by padding with a fill value (0 for black).
        neighbor_sum = convolve2d(img, kernel, mode='same', boundary='fill', fillvalue=0)
        
        # Apply the majority rule: pixels become 1 if the sum is over the threshold, else 0.
        img = (neighbor_sum > majority_threshold).astype(int)

    return img

# --- Simulation Parameters ---
# Image size
N = 20
# Neighborhood radius (delta << N)
delta = 1 
# Number of iterations (t)
num_iterations = 5

# --- Create Initial Image ---
# We create an "L" shape to show the evolution of both convex and concave corners.
initial_image = np.zeros((N, N), dtype=int)
# Vertical bar of the L
initial_image[5:15, 5:8] = 1
# Horizontal bar of the L
initial_image[12:15, 5:12] = 1

# --- Run Simulation and Print Result ---
final_image = run_median_filter_simulation(initial_image, delta, num_iterations)

print(f"\nFinal Image after {num_iterations} iterations:")
# Print the final image using the same clear format
for row in final_image:
    print("".join(map(str, row)).replace('0', '.').replace('1', '#'))