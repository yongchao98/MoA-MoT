import numpy as np
from scipy.ndimage import median_filter

def print_image_ascii(image):
    """Prints a binary numpy array to the console using ASCII characters."""
    for row in image:
        print(" ".join(['##' if pixel > 0 else '..' for pixel in row]))
    print("\n")

def run_median_filter_simulation():
    """
    Demonstrates the long-term effect of an iterative local median filter
    on a binary image.
    """
    # 1. Setup the image and filter parameters
    N = 40  # Image size: N x N
    delta = 2  # Radius of the local median filter
    iterations = 50 # Number of times to apply the filter

    # Create the initial image: a hollow square
    # This shape has an outer (convex) boundary and an inner (concave) boundary.
    image = np.zeros((N, N), dtype=int)
    image[10:30, 10:30] = 1  # White square
    image[14:26, 14:26] = 0  # Black hole in the middle

    # 2. Create the circular structuring element for the median filter
    # This defines the neighborhood ||x-y|| <= delta
    y, x = np.ogrid[-delta:delta+1, -delta:delta+1]
    footprint = x**2 + y**2 <= delta**2

    # 3. Print the initial state
    print("="*40)
    print("INITIAL IMAGE (t=0)")
    print("="*40)
    print_image_ascii(image)

    # 4. Apply the filter iteratively
    filtered_image = image.copy()
    for t in range(iterations):
        # Stop if the image becomes uniform (all black or all white)
        if np.all(filtered_image == 0) or np.all(filtered_image == 1):
            print(f"Image stabilized at iteration t={t}.\n")
            break
        filtered_image = median_filter(filtered_image, footprint=footprint, mode='constant', cval=0)

    # 5. Print the final state
    print("="*40)
    print(f"FINAL IMAGE (after {iterations} iterations or stabilization)")
    print("="*40)
    print_image_ascii(filtered_image)

    # 6. Explain the result
    print("="*40)
    print("EXPLANATION")
    print("="*40)
    print("The iterative application of a local median filter causes the edges of shapes to move based on their curvature.")
    print("\nKey observations:")
    print("1. Smoothing: Sharp corners on the square were rounded off.")
    print("2. Shrinking (Erosion): The outer, convex boundary of the white frame moved inwards.")
    print("3. Growth (Dilation): The inner, concave boundary of the white frame moved outwards, filling the hole.")
    print("\nThis process is an approximation of 'motion by mean curvature', which acts to minimize the total length of edges in the image.")
    print("As t -> infinity, this causes any finite shape to shrink and ultimately vanish, leading to a uniform image (in this case, all black).")

if __name__ == '__main__':
    run_median_filter_simulation()