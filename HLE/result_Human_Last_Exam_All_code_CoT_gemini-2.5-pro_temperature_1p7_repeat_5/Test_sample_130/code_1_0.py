import numpy as np
from scipy.ndimage import convolve

def apply_majority_filter(image, size=3):
    """
    Applies a majority filter to a binary image using a square kernel of a given size.
    The new pixel value is 1 if the number of 1s in the neighborhood is strictly
    greater than half the kernel size, and 0 otherwise.
    
    Args:
        image (np.array): The input binary image (0s and 1s).
        size (int): The odd integer size of the square neighborhood (e.g., 3 for a 3x3 kernel).
        
    Returns:
        np.array: The image after applying the filter.
    """
    if size % 2 == 0:
        raise ValueError("Kernel size must be odd.")
        
    # The convolve function sums the values in the neighborhood.
    kernel = np.ones((size, size))
    neighbor_sum = convolve(image.astype(float), kernel, mode='constant', cval=0)
    
    # The threshold for the majority is half the number of pixels in the kernel.
    threshold = (size * size) / 2.0
    
    # Return a new binary image based on the threshold.
    return (neighbor_sum > threshold).astype(image.dtype)

def run_simulation():
    """
    Runs a simulation to show the effect of repeated median filtering on a shape's edges.
    """
    # Create an 10x10 image with a 6x6 white square in the middle.
    image_size = 12
    shape_size = 6
    offset = (image_size - shape_size) // 2
    
    img = np.zeros((image_size, image_size), dtype=np.uint8)
    img[offset:offset+shape_size, offset:offset+shape_size] = 1
    
    print("--- Simulation of Repeated Local Median Filtering ---")
    print("The local median on a binary image is a majority filter.")
    print("We will apply a 3x3 filter repeatedly to a white square.\n")
    
    # Print the initial image
    print(f"Initial Image (t=0):")
    print(img)
    print("-" * 25)

    # Apply the filter for a few iterations and print the result each time.
    num_iterations = 6
    for t in range(1, num_iterations + 1):
        img = apply_majority_filter(img, size=3)
        area = np.sum(img)
        print(f"Image after t={t} iteration(s):")
        print(img)
        if area == 0:
            print("\nThe shape has vanished completely.")
            break
        else:
            print("-" * 25)

if __name__ == '__main__':
    run_simulation()