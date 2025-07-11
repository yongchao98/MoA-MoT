import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter

# This script requires the scikit-image library to draw the initial shape.
# You can install it via pip:
# pip install scikit-image

def create_star_image(size=200, points=5):
    """Creates a binary image with a star shape."""
    try:
        import skimage.draw
    except ImportError:
        print("Scikit-image not found. Please install it using 'pip install scikit-image'")
        print("Creating a simple square image instead.")
        image = np.zeros((size, size), dtype=np.uint8)
        image[size//4:-size//4, size//4:-size//4] = 1
        return image

    image = np.zeros((size, size), dtype=np.uint8)
    radius_outer = size / 2.5
    radius_inner = radius_outer / 2.5
    center_y, center_x = size // 2, size // 2
    
    # Generate vertices for the star polygon
    angles = np.linspace(0, 2 * np.pi, 2 * points + 1)
    radii = np.tile([radius_outer, radius_inner], points)
    
    poly_y = center_y + radii * np.sin(angles[:-1])
    poly_x = center_x + radii * np.cos(angles[:-1])
    
    # Draw the polygon on the image
    rr, cc = skimage.draw.polygon(poly_y, poly_x)
    rr = np.clip(rr, 0, size - 1)
    cc = np.clip(cc, 0, size - 1)
    image[rr.astype(int), cc.astype(int)] = 1
    return image

def main():
    """
    Main function to run the simulation and display the results.
    """
    # --- Parameters ---
    image_size = 200        # N: Height and width of the image
    neighborhood_radius = 3 # δ: Radius of the neighborhood for the median filter
    iterations_to_show = [0, 5, 20, 100] # Time steps t to visualize

    # The filter size is the width of the square neighborhood
    filter_size = 2 * neighborhood_radius + 1

    print("--- Simulating Edge Evolution with a Repeated Median Filter ---")
    print(f"Image Size (N x N): {image_size} x {image_size}")
    print(f"Neighborhood Radius (δ): {neighborhood_radius}")
    print(f"Median Filter Size: {filter_size} x {filter_size}")
    print("\nThe simulation shows what happens to the edges of a star shape over time.")
    print("Observe how sharp corners are rounded and the entire shape shrinks.")

    # Create the initial binary image
    image = create_star_image(image_size)
    
    # Set up the plot to show evolution over time
    num_plots = len(iterations_to_show)
    fig, axes = plt.subplots(1, num_plots, figsize=(num_plots * 4, 4.5))
    fig.suptitle("Evolution of Edges Under Repeated Local Median Filtering", fontsize=16)

    current_image = image.copy()
    plot_index = 0
    max_iter = max(iterations_to_show)

    # Main loop: Apply the filter repeatedly
    for t in range(max_iter + 1):
        if t in iterations_to_show:
            ax = axes[plot_index]
            ax.imshow(current_image, cmap='gray', vmin=0, vmax=1)
            ax.set_title(f"t = {t} iterations")
            ax.set_xticks([])
            ax.set_yticks([])
            plot_index += 1
        
        # Apply the local median filter
        current_image = median_filter(current_image, size=filter_size)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

    print("\nConclusion: As t -> ∞, the edges of a closed shape are smoothed and eventually eliminated.")

if __name__ == '__main__':
    main()
