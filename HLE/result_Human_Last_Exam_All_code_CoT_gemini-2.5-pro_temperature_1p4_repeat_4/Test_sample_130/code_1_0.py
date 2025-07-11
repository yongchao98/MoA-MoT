import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter

def simulate_median_filter_evolution():
    """
    Simulates and visualizes the evolution of a binary image
    under iterative application of a local median filter.
    """
    # --- Parameters ---
    N = 100  # Image size: N x N
    delta = 2  # Radius of the median filter neighborhood
    filter_size = 2 * delta + 1
    iterations_to_show = [0, 5, 20, 100]

    # --- Print the "Equation" and parameters as requested ---
    print("This simulation demonstrates the iterative application of the local median filter.")
    print("The evolution rule is defined by the equation: f_{t+1}(x) = Median_{||y-x|| <= delta}(f_t(y))")
    print("\nParameters used in the simulation:")
    print(f"Image Size (N): {N}")
    print(f"Neighborhood Radius (delta): {delta}")
    print(f"Resulting Filter Size ({2*delta+1}x{2*delta+1}): {filter_size}")
    print("-" * 30)
    print("We will start with a cross shape and observe how its edges evolve.")

    # --- Create the initial image (a cross shape) ---
    image = np.zeros((N, N), dtype=int)
    cross_width = 8
    center = N // 2
    arm_length = N // 3
    # Horizontal arm
    image[center - cross_width // 2 : center + cross_width // 2, center - arm_length : center + arm_length] = 1
    # Vertical arm
    image[center - arm_length : center + arm_length, center - cross_width // 2 : center + cross_width // 2] = 1

    # --- Setup visualization ---
    num_plots = len(iterations_to_show)
    fig, axes = plt.subplots(1, num_plots, figsize=(num_plots * 4, 4.5))
    fig.suptitle("Evolution of a Shape under Iterative Median Filtering", fontsize=16)

    current_image = image.copy()
    plot_idx = 0

    # --- Run the simulation and plot at specified intervals ---
    for i in range(max(iterations_to_show) + 1):
        # Check if the current iteration is one we want to display
        if i in iterations_to_show:
            ax = axes[plot_idx]
            ax.imshow(current_image, cmap='gray', interpolation='none')
            ax.set_title(f"t = {i} iterations")
            ax.set_xticks([])
            ax.set_yticks([])
            plot_idx += 1
        
        # Apply the median filter for the next iteration
        current_image = median_filter(current_image, size=filter_size)

    # --- Print Final Explanation ---
    print("\nObservation from the plots:")
    print("As the number of iterations 't' increases, the sharp concave and convex corners of the cross are smoothed out.")
    print("The shape evolves to reduce its boundary length, tending towards a circle.")
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

# Run the simulation
simulate_median_filter_evolution()
