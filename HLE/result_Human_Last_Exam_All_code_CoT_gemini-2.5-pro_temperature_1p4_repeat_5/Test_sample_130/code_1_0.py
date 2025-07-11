import numpy as np
from scipy.ndimage import generic_filter

def create_disk_footprint(radius):
    """Creates a circular structuring element (footprint) for neighborhood operations."""
    y, x = np.ogrid[-radius:radius+1, -radius:radius+1]
    # Create a mask for a disk shape
    mask = x**2 + y**2 <= radius**2
    return mask

def print_image(image, title):
    """Helper function to print the image grid."""
    print(title)
    # The problem requests printing each number in the final equation,
    # which we interpret as printing the final state of the image matrix.
    for row in image:
        # We replace 0 with '.' and 1 with '#' for better visibility
        printable_row = "".join(['#' if p == 1 else '.' for p in row])
        print(printable_row)
    print("-" * (len(image[0]) + 1))

def main():
    """
    Demonstrates the effect of an iterative local median filter on a binary image.
    """
    # Define parameters
    N = 25  # Image size
    delta = 2  # Radius for the median filter
    iterations_to_show = [0, 1, 3, 6, 10]

    # Create an initial binary image with a non-convex 'L' shape
    image = np.zeros((N, N), dtype=int)
    image[6:19, 7:11] = 1
    image[15:19, 11:17] = 1

    # Create the disk-shaped footprint for the filter
    footprint = create_disk_footprint(delta)

    # Apply the filter iteratively
    current_image = image.copy()
    for t in range(max(iterations_to_show) + 1):
        if t in iterations_to_show:
            print_image(current_image, f"Image after t={t} iterations:")

        # The core operation: applying the local median filter
        # It calculates the median in a disk neighborhood for each pixel.
        current_image = generic_filter(current_image, np.median, footprint=footprint, mode='constant', cval=0)

if __name__ == '__main__':
    main()
