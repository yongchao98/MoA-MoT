import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter
from skimage.morphology import disk
import warnings

# Suppress a specific warning from matplotlib that might occur in some environments
# when creating many figures, not critical for the logic.
warnings.filterwarnings("ignore", message="More than 20 figures have been opened.")

def run_simulation_and_explain():
    """
    This function runs a simulation to demonstrate the effect of an iterative
    local median filter on a binary image and prints a detailed explanation.
    """
    # 1. Define image and filter parameters
    N = 100  # Image size: N x N
    delta = 3  # Radius of the local median filter

    # 2. Create an initial binary image with a cross shape
    # This shape has both convex and concave corners, making it a good example.
    image = np.zeros((N, N), dtype=np.uint8)
    image[N//2 - 25 : N//2 + 25, N//2 - 8 : N//2 + 8] = 1
    image[N//2 - 8 : N//2 + 8, N//2 - 25 : N//2 + 25] = 1

    # 3. Define the structuring element for the median filter (a disk)
    # This corresponds to the ||x-y|| <= delta neighborhood.
    footprint = disk(delta)

    # 4. Set up the simulation for iterative application
    # We will show the image at these specific iteration counts.
    iterations_to_show = [0, 10, 30, 70]
    images_over_time = []
    current_image = image.copy()

    max_iter = max(iterations_to_show)
    for t in range(max_iter + 1):
        if t in iterations_to_show:
            images_over_time.append((t, current_image.copy()))
        
        # Apply the local median filter
        current_image = median_filter(current_image, footprint=footprint)

    # 5. Display the visual results of the simulation
    fig, axes = plt.subplots(1, len(images_over_time), figsize=(16, 5))
    if len(images_over_time) == 1: # Handle case of a single subplot
        axes = [axes]
    for ax, (t, img) in zip(axes, images_over_time):
        ax.imshow(img, cmap='gray', vmin=0, vmax=1)
        ax.set_title(f'Iteration t = {t}')
        ax.axis('off')

    fig.suptitle(f'Evolution of Edges Under Iterative Median Filtering (N={N}, δ={delta})', fontsize=16)
    plt.show()

    # 6. Print the detailed explanation
    print("--- Analysis of Iterative Median Filtering on Binary Image Edges ---")
    print("\nThe simulation above shows a binary image at different stages of a process where a local median filter is applied repeatedly.")
    print(f"\nParameters Used: Image Size = {N}x{N}, Filter Radius (δ) = {delta}\n")
    print("The question is: What happens to the edges of the image as t -> infinity?\n")
    print("1. SMOOTHING OF CORNERS:")
    print("   The first thing to happen is that sharp corners (areas of high curvature) are smoothed out.")
    print("   A pixel at a sharp convex corner is surrounded by a majority of pixels from the background, so it flips color, eroding the corner.")
    print("   Conversely, a sharp concave corner gets filled in.")
    
    print("\n2. MOTION BY MEAN CURVATURE:")
    print("   This process causes the boundary of a shape to move inward, perpendicular to the edge itself.")
    print("   The speed of this movement is proportional to the local curvature of the edge.")
    print("   - Convex parts (like the outside of the cross) shrink.")
    print("   - Concave parts (like the inside corners of the cross) are filled in.")
    print("   This causes the shape to become more rounded and circular over time.")
    
    print("\n3. SHRINKING AND DISAPPEARANCE:")
    print("   As t -> infinity, this shrinking process continues. Any isolated shape with a closed boundary will continuously shrink.")
    print("   Eventually, the shape will shrink to a single point and then vanish completely.")
    
    print("\nCONCLUSION:")
    print("   As t -> infinity, the edges of any closed shape will smooth out, shrink, and ultimately disappear.")
    print("   The final state of the image will be monochromatic (either all 0s or all 1s), as all features are eroded away.")

# Run the simulation and print the explanation
if __name__ == '__main__':
    run_simulation_and_explain()