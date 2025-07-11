import numpy as np
from scipy.ndimage import median_filter

def simulate_edge_evolution(N=20, delta=1, iterations=6):
    """
    Simulates the effect of an iterative local median filter on a binary image.

    Args:
        N (int): The size of the N x N image.
        delta (int): The radius of the median filter neighborhood.
        iterations (int): The number of times to apply the filter.
    """
    # Create a binary image with a white square on a black background
    image = np.zeros((N, N), dtype=int)
    start, end = N // 4, 3 * N // 4
    image[start:end, start:end] = 1

    print("This script demonstrates what happens to the edges of a shape in a binary image")
    print("when a local median filter is applied repeatedly.")
    print("We start with a white square on a black background.")
    print("1 represents a white pixel, and 0 represents a black pixel.\n")

    print(f"Initial Image (t=0):")
    for row in image:
        print(' '.join(map(str, row)))
    print("-" * (2 * N))

    current_image = image.copy()
    # The filter size is (2*delta + 1) x (2*delta + 1)
    filter_size = 2 * delta + 1

    for t in range(1, iterations + 1):
        # Apply the median filter
        # The 'mode' and 'cval' handle pixels at the image border.
        # 'constant' with 'cval=0' means we assume the image is surrounded by black pixels.
        current_image = median_filter(current_image, size=filter_size, mode='constant', cval=0)

        print(f"Image after t={t} iteration(s):")
        # Print the image array to show the changes
        for row in current_image:
            print(' '.join(map(str, row)))

        # Check if the shape has disappeared
        if np.sum(current_image) == 0:
            print("\nThe shape has completely disappeared.")
            break
        print("-" * (2 * N))


if __name__ == '__main__':
    # Run the simulation. A small delta (like 1) is sufficient to see the effect.
    simulate_edge_evolution(N=20, delta=1, iterations=10)