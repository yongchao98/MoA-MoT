import numpy as np
from scipy.ndimage import median_filter

def simulate_median_filtering(image_size, rect_shape, delta, iterations):
    """
    Simulates the effect of a repeated local median filter on a binary image.

    Args:
        image_size (int): The width and height of the square image (N).
        rect_shape (tuple): The (height, width) of the initial white rectangle.
        delta (int): The radius for the median filter (δ).
        iterations (int): The number of times to apply the filter (simulating t).
    """
    # 1. Create a binary black and white image with a rectangle in the center
    image = np.zeros((image_size, image_size), dtype=int)
    rect_h, rect_w = rect_shape
    start_row = (image_size - rect_h) // 2
    start_col = (image_size - rect_w) // 2
    image[start_row : start_row + rect_h, start_col : start_col + rect_w] = 1

    print("--- Initial Image ---")
    print("A '1' represents a white pixel, and a '0' represents a black pixel.")
    print(image)
    print("\n")

    # The size of the median filter neighborhood corresponding to radius delta
    # For a square neighborhood, the side length is 2*delta + 1
    filter_size = 2 * delta + 1
    
    # 2. Apply the local median filter repeatedly
    filtered_image = image.copy()
    for i in range(iterations):
        filtered_image = median_filter(filtered_image, size=filter_size, mode='constant', cval=0)

    print(f"--- Image after {iterations} iterations with delta={delta} ---")
    print("Notice the sharp corners have been rounded.")
    print("The final matrix of pixel values is:")
    print(filtered_image)


if __name__ == '__main__':
    # Parameters for the simulation
    N = 40      # Image size
    DELTA = 4   # Radius of the median operator
    T = 20      # Number of iterations to approximate t -> infinity
    
    # Shape of the initial white object
    rectangle_shape = (20, 20)

    # Check constraints
    if DELTA >= N:
        print("Error: Delta (δ) should be much smaller than N.")
    else:
        simulate_median_filtering(N, rectangle_shape, DELTA, T)
