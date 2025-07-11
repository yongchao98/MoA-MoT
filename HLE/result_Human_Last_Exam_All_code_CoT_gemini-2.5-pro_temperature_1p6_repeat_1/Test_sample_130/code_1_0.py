import numpy as np
from scipy.ndimage import median_filter

def print_image(image, title):
    """Helper function to print the image array with a title."""
    print(title)
    # Replace 0s with '.' and 1s with '#' for better visualization
    for row in image:
        print(" ".join(['#' if pixel == 1 else '.' for pixel in row]))
    print("-" * (2 * image.shape[1]))

def demonstrate_median_filter_evolution():
    """
    Demonstrates the effect of a repeated local median filter on a binary image.
    """
    # Parameters
    image_size = 12
    # The radius of the filter is 1, meaning a 3x3 neighborhood
    filter_size = 3
    iterations = 5

    # 1. Create a binary image with an L-shape
    # This shape has convex (outer) corners and one concave (inner) corner
    image = np.zeros((image_size, image_size), dtype=int)
    image[2:9, 2:5] = 1  # Vertical part of 'L'
    image[7:9, 5:9] = 1  # Horizontal part of 'L'

    # 2. Print the initial image
    print_image(image, "Initial Image (t=0):")

    # 3. Apply the median filter repeatedly
    evolved_image = image.copy()
    for _ in range(iterations):
        evolved_image = median_filter(evolved_image, size=filter_size, mode='constant', cval=0)
    
    # 4. Print the final image
    print(f"\nAfter {iterations} iterations, the sharp corners are smoothed.\n"
          "The convex (outer) corners have been eroded, and the\n"
          "concave (inner) corner has started to fill in.\n")
    print_image(evolved_image, f"Image after {iterations} iterations:")

if __name__ == '__main__':
    demonstrate_median_filter_evolution()