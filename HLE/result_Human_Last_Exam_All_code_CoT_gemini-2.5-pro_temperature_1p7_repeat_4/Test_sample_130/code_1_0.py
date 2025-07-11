import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter

def plot_image_evolution(initial_image, iterations, filter_size):
    """
    Applies a median filter iteratively and plots the evolution of the image.

    Args:
        initial_image (np.array): The starting binary image.
        iterations (list): A list of iteration numbers to plot.
        filter_size (int): The size of the median filter neighborhood.
    """
    img = initial_image.copy()
    
    plot_stages = {0: initial_image}
    
    print(f"Applying a {filter_size}x{filter_size} median filter iteratively...")

    for i in range(1, max(iterations) + 1):
        img = median_filter(img, size=filter_size)
        if i in iterations:
            plot_stages[i] = img.copy()
            print(f"Captured image state at iteration t = {i}")

    fig, axes = plt.subplots(1, len(plot_stages), figsize=(15, 5))
    sorted_stages = sorted(plot_stages.items())

    for ax, (i, stage_img) in zip(axes, sorted_stages):
        ax.imshow(stage_img, cmap='gray', vmin=0, vmax=1)
        ax.set_title(f'Iteration t = {i}')
        ax.set_xticks([])
        ax.set_yticks([])
    
    plt.suptitle("Evolution of Edges Under Iterative Median Filtering", fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

# --- Main Execution ---
if __name__ == "__main__":
    # Define image and shape parameters
    N = 100  # Image size (N x N)
    delta = 7 # Corresponds to filter size, must be odd
    
    # Create an initial binary image with a cross shape
    # This shape has many sharp convex and concave corners
    image = np.zeros((N, N), dtype=int)
    cross_width = 8
    center = N // 2
    arm_length = N // 3
    
    image[center - cross_width // 2 : center + cross_width // 2, center - arm_length : center + arm_length] = 1
    image[center - arm_length : center + arm_length, center - cross_width // 2 : center + cross_width // 2] = 1

    # Define the iteration numbers at which to view the image
    evolution_iterations = [1, 5, 15, 30]

    # Run the simulation and plot the results
    plot_image_evolution(image, evolution_iterations, filter_size=delta)