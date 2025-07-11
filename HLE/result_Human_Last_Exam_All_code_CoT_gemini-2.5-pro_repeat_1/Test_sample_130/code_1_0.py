import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter

def run_median_filter_simulation():
    """
    Simulates the long-term application of a local median filter on a binary image.
    """
    # 1. Define parameters
    image_size = 200  # N
    delta = 5         # Radius for the median filter neighborhood
    iterations = [0, 10, 50, 100] # Time steps t to visualize

    # The problem defines the local median function as f_{t+delta}(x) = Median_{||x-y||<=delta}
    # We will print the parameters for each step of the simulation.
    print(f"Simulation Parameters:")
    print(f"Image Size (N x N): {image_size}x{image_size}")
    print(f"Neighborhood Radius (delta): {delta}")
    print("-" * 20)

    # 2. Create the initial binary image (a donut shape)
    # This shape has a convex outer boundary and a concave inner boundary.
    image = np.zeros((image_size, image_size), dtype=np.uint8)
    y, x = np.ogrid[-image_size/2:image_size/2, -image_size/2:image_size/2]
    
    # The distance of each pixel from the center
    dist_from_center = np.sqrt(x**2 + y**2)

    # Radii for the donut shape
    outer_radius = 70
    inner_radius = 40
    
    # Create the donut mask
    donut_mask = (dist_from_center >= inner_radius) & (dist_from_center <= outer_radius)
    image[donut_mask] = 1

    # 3. Create the structuring element (footprint) for the median filter
    # This defines the neighborhood ||x-y|| <= delta
    y_foot, x_foot = np.ogrid[-delta:delta+1, -delta:delta+1]
    footprint = x_foot**2 + y_foot**2 <= delta**2

    # 4. Run the simulation and visualize
    fig, axes = plt.subplots(1, len(iterations), figsize=(16, 5))
    fig.suptitle(f'Evolution of Edges under Iterative Median Filtering (delta={delta})', fontsize=16)

    current_image = image.copy()
    
    # The problem specifies the operator f_{t+delta}, implying time steps.
    # We simulate this by iteratively applying the filter.
    simulation_step = 0
    for i, t in enumerate(range(max(iterations) + 1)):
        if t in iterations:
            ax = axes[simulation_step]
            ax.imshow(current_image, cmap='gray', vmin=0, vmax=1)
            ax.set_title(f'After t = {t} iterations')
            ax.set_xticks([])
            ax.set_yticks([])
            print(f"Displaying image at t = {t}")
            simulation_step += 1
        
        # Apply the local median filter
        current_image = median_filter(current_image, footprint=footprint)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

if __name__ == '__main__':
    run_median_filter_simulation()