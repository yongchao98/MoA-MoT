# Note: This code requires the 'scipy' library.
# You can install it using: pip install scipy
import numpy as np
from scipy.ndimage import median_filter

def demonstrate_edge_evolution():
    """
    Demonstrates the effect of a repeated local median filter on the edges
    of a binary image.
    """
    # The local median function is defined as: f_{t+δ}(x) = Median_{||x-y||<=δ}
    # We will simulate its application over several steps of t.

    # 1. Define the parameters for the function.
    delta = 1  # δ: radius of the neighborhood
    filter_size = 2 * delta + 1 # The filter neighborhood is a square of this size.

    # 2. Create an initial binary image (I at t=0) with a non-convex shape.
    #    Here, 0 is black and 1 is white. The shape is a black 'U'.
    image_t = np.array([
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1],
        [1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1],
        [1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1],
        [1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1],
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    ], dtype=int)

    num_iterations = 5

    print(f"Demonstrating the operator f with delta = {delta} (a {filter_size}x{filter_size} median filter).")
    print("We observe the image state at each time step t.\n")

    print(f"Image at t = 0 (Initial state):\n{image_t}\n" + "-"*40)

    # 3. Iterate from t=1 to num_iterations, applying the filter each time.
    for t in range(1, num_iterations + 1):
        # Apply the median filter. mode='constant' defines the pixel value
        # outside the boundary, which is set to white (cval=1).
        image_t = median_filter(image_t, size=filter_size, mode='constant', cval=1)

        # Output the numbers in the image at the current step t.
        print(f"Image at t = {t}:\n{image_t}\n" + "-"*40)

    print("Observe how the sharp outer corners (convex) become rounded and the")
    print("inner part of the 'U' (concave) fills in, smoothing the boundary.")

demonstrate_edge_evolution()