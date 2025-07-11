import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import generic_filter

def simulate_edge_evolution():
    """
    Simulates the effect of a repeated local median filter on a binary image.
    """
    # 1. Define parameters
    N = 100  # Image size: N x N
    delta = 3  # Radius of the local median filter's neighborhood
    iterations = [0, 5, 20, 60] # Timesteps to visualize

    # 2. Create the initial binary image (a cross shape)
    # This shape has both convex and concave corners to show the effect clearly.
    image = np.zeros((N, N), dtype=int)
    image[N//2 - 20 : N//2 + 20, N//2 - 5 : N//2 + 5] = 1
    image[N//2 - 5 : N//2 + 5, N//2 - 20 : N//2 + 20] = 1

    # 3. Define the local median filter (as a majority filter)
    # The filter function returns 1 if the majority of pixels in the
    # neighborhood are 1, and 0 otherwise. This is the essence of
    # a median filter on binary data.
    def majority_rule_filter(neighborhood):
        # Tie-breaking: if mean is exactly 0.5, it becomes 0.
        return 1 if np.mean(neighborhood) > 0.5 else 0

    # Create a circular footprint for the neighborhood
    y, x = np.ogrid[-delta:delta+1, -delta:delta+1]
    footprint = x**2 + y**2 <= delta**2

    # 4. Run the simulation and visualize the results
    fig, axes = plt.subplots(1, len(iterations), figsize=(16, 5))
    fig.suptitle(f'Edge Evolution under Repeated Local Median Filter (delta={delta})', fontsize=16)

    current_image = image.copy()
    iter_idx = 0

    for t in range(max(iterations) + 1):
        if t in iterations:
            ax = axes[iter_idx]
            ax.imshow(current_image, cmap='gray', interpolation='nearest')
            ax.set_title(f't = {t}')
            ax.set_xticks([])
            ax.set_yticks([])
            iter_idx += 1

        # Apply the filter to get the next state of the image
        current_image = generic_filter(current_image, majority_rule_filter, footprint=footprint, mode='constant', cval=0)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    print("Displaying plot of image evolution. Close the plot window to exit.")
    plt.show()

# Execute the simulation
simulate_edge_evolution()