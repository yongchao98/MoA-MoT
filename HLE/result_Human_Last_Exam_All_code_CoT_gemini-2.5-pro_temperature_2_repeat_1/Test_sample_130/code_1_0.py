import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter

def create_cross_image(N=100, arm_thickness=10, arm_length=30):
    """Creates a binary image with a cross shape."""
    image = np.zeros((N, N), dtype=int)
    center = N // 2
    
    # Horizontal arm
    image[center - arm_thickness//2 : center + arm_thickness//2, 
          center - arm_length//2 : center + arm_length//2] = 1
    # Vertical arm
    image[center - arm_length//2 : center + arm_length//2,
          center - arm_thickness//2 : center + arm_thickness//2] = 1
          
    return image

def run_simulation():
    """
    Runs a simulation of applying a median filter iteratively to a binary image.
    """
    # Parameters
    N = 100 # Image size
    delta = 3  # Radius of the median filter neighborhood
    filter_size = 2 * delta + 1
    iterations = 30
    
    # Create the initial image
    image = create_cross_image(N=N, arm_thickness=20, arm_length=70)
    initial_white_pixels = np.sum(image)
    
    print(f"Initial image size: {N}x{N}")
    print(f"Median filter neighborhood radius δ: {delta} (size: {filter_size}x{filter_size})")
    print(f"Running simulation for {iterations} iterations.\n")
    
    # Set up plots
    # We will show the image at a few key steps
    display_steps = [0, 5, 10, 15, 20, 25]
    fig, axes = plt.subplots(1, len(display_steps), figsize=(18, 4))
    
    # Run the simulation
    current_image = image.copy()
    
    for i in range(iterations + 1):
        if i in display_steps:
            ax_index = display_steps.index(i)
            axes[ax_index].imshow(current_image, cmap='gray', vmin=0, vmax=1)
            axes[ax_index].set_title(f"t = {i}")
            axes[ax_index].axis('off')

        white_pixels = np.sum(current_image)
        # We must output the numbers as requested.
        print(f"Iteration t={i:2d}: Number of white pixels = {white_pixels:4d} ({(white_pixels/initial_white_pixels)*100:.1f}%)")

        if white_pixels == 0:
            print("\nImage has become all black. Halting simulation.")
            # Fill remaining plots if necessary
            for j in range(ax_index + 1, len(display_steps)):
                axes[j].imshow(current_image, cmap='gray', vmin=0, vmax=1)
                axes[j].set_title(f"t >= {i}")
                axes[j].axis('off')
            break
            
        # Apply the local median filter
        current_image = median_filter(current_image, size=filter_size)
        
    fig.suptitle("Evolution of a Binary Image under Iterative Median Filtering", fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

# Execute the simulation
run_simulation()
