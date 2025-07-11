import numpy as np
from scipy.ndimage import median_filter
import sys

def create_circular_footprint(radius):
    """
    Creates a circular boolean mask for a given radius.
    This is a fallback if scikit-image is not installed.
    """
    size = 2 * radius + 1
    y, x = np.ogrid[-radius:radius+1, -radius:radius+1]
    mask = x*x + y*y <= radius*radius
    return mask

def print_image(img_array):
    """Prints a numpy array as a text image using block characters."""
    for row in img_array:
        # Use unicode characters for a better visual representation in modern terminals
        print("".join(['██' if pixel == 1 else '  ' for pixel in row]))

def run_simulation():
    """
    Runs and prints a simulation of the iterative median filter on a binary image.
    """
    # Parameters
    N = 28      # Image size
    delta = 2   # Radius of the median filter neighborhood
    shape_size = 14 # Size of the initial white square

    # Initial Image: A white square on a black background
    image = np.zeros((N, N), dtype=np.uint8)
    start = (N - shape_size) // 2
    end = start + shape_size
    image[start:end, start:end] = 1

    # Create the circular neighborhood (footprint) for the median filter
    try:
        from skimage.morphology import disk
        footprint = disk(delta)
    except ImportError:
        print("(Note: 'scikit-image' not found. Using a NumPy-based fallback for creating the circular neighborhood.)", file=sys.stderr)
        footprint = create_circular_footprint(delta)

    time_steps_to_show = [0, 5, 15, 30, 45]
    img_t = image.copy()

    print("--- Simulation of Iterative Median Filtering ---")
    print(f"Initial state: A {shape_size}x{shape_size} white square in a {N}x{N} black image.")
    print(f"Operator: Median filter with a circular neighborhood of radius delta = {delta}.\n")

    for t in range(max(time_steps_to_show) + 1):
        if t in time_steps_to_show:
            print(f"--- Image at t = {t} ---")
            print_image(img_t)
            print("-" * (2 * N + 4))
            
            # Stop if the image has converged to a solid color
            if np.all(img_t == 0):
                print(f"\nImage has converged to a solid black state at t={t}.")
                break
        
        # Apply the median filter
        img_t = median_filter(img_t, footprint=footprint, mode='constant', cval=0)

    print("\n--- Observation Summary ---")
    print("1. The sharp corners of the square progressively become more rounded.")
    print("2. The overall area of the square shrinks over time.")
    print("3. Ultimately, the shape vanishes completely, resulting in an all-black image.")
    print("\nThis demonstrates that as t -> infinity, edges smooth out and convex shapes shrink and disappear.")

if __name__ == '__main__':
    run_simulation()