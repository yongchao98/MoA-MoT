import numpy as np
from scipy.ndimage import median_filter

def demonstrate_iterative_median_filter():
    """
    Demonstrates the effect of an iterative median filter on a binary image.
    This simulation shows how edges evolve over time.
    """
    # 1. Create an initial binary image with a shape (a cross)
    N = 15
    image = np.zeros((N, N), dtype=int)
    center = N // 2
    cross_width = 0 # Use a width of 1 pixel (indices are inclusive)
    cross_len = 5
    # Create vertical bar of the cross
    image[center - cross_len:center + cross_len + 1, center - cross_width:center + cross_width + 1] = 1
    # Create horizontal bar of the cross
    image[center - cross_width:center + cross_width + 1, center - cross_len:center + cross_len + 1] = 1

    # 2. Define filter parameters
    delta = 1
    filter_size = 2 * delta + 1
    iterations_to_show = [0, 1, 3, 6, 9]

    print("This script simulates the iterative application of a local median filter.")
    print(f"The initial image is a {N}x{N} grid with a cross shape.")
    print(f"The median filter has a radius delta = {delta} (a {filter_size}x{filter_size} neighborhood).")
    print("-" * (2 * N))

    current_image = image.copy()

    # 3. Iterate and print the image at different steps
    for i in range(max(iterations_to_show) + 1):
        if i in iterations_to_show:
            print(f"Image after t = {i} iterations:")
            # Note: There is no equation in this problem, so the instruction to
            # "output each number in the final equation" is not applicable here.
            for row in current_image:
                # Print a more visual representation of the image
                print(" ".join(['■' if pixel == 1 else '·' for pixel in row]))
            print("-" * (2 * N))
            
            # Stop if the image has no more shapes
            if np.sum(current_image) == 0:
                print(f"The shape has completely disappeared after {i} iterations.")
                break

        # Apply the median filter for the next iteration
        current_image = median_filter(current_image, size=filter_size, mode='constant', cval=0)

if __name__ == '__main__':
    demonstrate_iterative_median_filter()