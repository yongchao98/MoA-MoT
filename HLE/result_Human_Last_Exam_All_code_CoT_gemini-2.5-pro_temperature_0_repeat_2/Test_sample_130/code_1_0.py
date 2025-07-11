import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter

def create_star_image(size=200):
    """Creates a binary image with a star shape."""
    image = np.zeros((size, size), dtype=np.uint8)
    center = size // 2
    
    # Create a polygon for the star
    from skimage.draw import polygon
    outer_radius = size / 2.5
    inner_radius = size / 5.0
    angles = np.linspace(0, 2 * np.pi, 11)
    outer_x = center + outer_radius * np.cos(angles)
    outer_y = center + outer_radius * np.sin(angles)
    inner_x = center + inner_radius * np.cos(angles + np.pi / 10)
    inner_y = center + inner_radius * np.sin(angles + np.pi / 10)
    
    poly_r = np.zeros(20)
    poly_c = np.zeros(20)
    
    poly_r[::2] = outer_y[:-1]
    poly_r[1::2] = inner_y[:-1]
    poly_c[::2] = outer_x[:-1]
    poly_c[1::2] = inner_x[:-1]
    
    rr, cc = polygon(poly_r, poly_c)
    image[rr, cc] = 1
    return image

def main():
    """
    Main function to demonstrate the effect of iterative median filtering.
    """
    # --- Parameters ---
    # The radius of the median filter neighborhood.
    delta = 5 
    # The number of times to apply the filter.
    num_iterations = 100

    # --- Setup ---
    # 1. Create the initial binary image with a complex shape.
    initial_image = create_star_image(size=200)
    
    # 2. Create the circular footprint for the median filter.
    # This defines the neighborhood ||x-y|| <= delta.
    y, x = np.ogrid[-delta:delta+1, -delta:delta+1]
    disk_footprint = x**2 + y**2 <= delta**2

    # --- Simulation ---
    # 3. Apply the median filter iteratively.
    print(f"Applying a median filter with radius delta = {delta} for {num_iterations} iterations...")
    
    # Start with the initial image
    processed_image = initial_image.copy()
    
    for i in range(num_iterations):
        processed_image = median_filter(processed_image, footprint=disk_footprint)
        # Optional: print progress
        if (i + 1) % 10 == 0:
            print(f"  ... iteration {i+1}/{num_iterations}")

    print("Simulation complete.")

    # --- Visualization ---
    # 4. Display the results.
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    
    axes[0].imshow(initial_image, cmap='gray')
    axes[0].set_title('Initial Image (t=0)')
    axes[0].axis('off')
    
    axes[1].imshow(processed_image, cmap='gray')
    axes[1].set_title(f'Final Image (t={num_iterations})')
    axes[1].axis('off')
    
    fig.suptitle(f'Effect of Iterative Median Filtering (delta={delta})', fontsize=16)
    plt.show()

if __name__ == '__main__':
    main()