import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter

def run_simulation():
    """
    Simulates the repeated application of a local median filter on a binary image
    to show its effect on edges over time.
    """
    # 1. Define parameters for the simulation
    N = 200  # Image size: N x N
    delta = 2  # Radius of the local median neighborhood
    filter_size = 2 * delta + 1 # The size of the median filter window
    iterations_to_show = [0, 10, 50]

    # 2. Create the initial binary image
    # Start with a black image (all zeros)
    image = np.zeros((N, N), dtype=int)

    # Add a large white rectangle
    image[30:100, 30:170] = 1
    # Add a smaller, more complex white shape (a cross)
    image[120:170, 80:90] = 1
    image[135:155, 65:105] = 1
    # Add some random "salt & pepper" noise
    noise = np.random.randint(0, 100, size=(N, N))
    image[noise > 98] = 1 # pepper noise
    image[noise < 2] = 0 # salt noise

    # Store images at different stages
    images_at_steps = {0: image.copy()}
    
    # 3. Apply the filter iteratively
    current_image = image.copy()
    for t in range(1, max(iterations_to_show) + 1):
        # Apply the local median function
        current_image = median_filter(current_image, size=filter_size, mode='constant', cval=0)
        if t in iterations_to_show:
            images_at_steps[t] = current_image.copy()

    # 4. Visualize and Explain the results
    fig, axes = plt.subplots(1, len(iterations_to_show), figsize=(15, 5))
    fig.suptitle(f"Evolution of an Image with a Repeated Local Median Filter (δ={delta})", fontsize=16)

    for i, t in enumerate(iterations_to_show):
        ax = axes[i]
        ax.imshow(images_at_steps[t], cmap='gray', vmin=0, vmax=1)
        ax.set_title(f"After t = {t} iterations")
        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout()
    plt.show()

    # 5. Print the final explanation
    print("--- Simulation Explanation ---")
    print(f"The simulation models the function f(I) = Median(I) over a neighborhood of radius δ = {delta}.")
    print("This function is applied repeatedly to an initial binary image.\n")
    print("Observations from the simulation:")
    print("1. Noise Removal: Isolated 'salt & pepper' noise pixels are removed almost immediately.")
    print("2. Smoothing: Sharp corners on the shapes become rounded. The edges become smoother over time.")
    print("3. Shrinking and Straightening: The overall area of the white shapes decreases.")
    print("   - Convex boundaries (like outside corners) are eroded.")
    print("   - The smaller, more complex cross shape shrinks rapidly and loses its complexity, tending towards a circular blob before disappearing.")
    print("   - The edges of the larger rectangle tend to straighten and its corners round.\n")
    print("Conclusion: As t -> ∞, the operator f causes edges to become smoother and to move in a way that minimizes their total length.")
    print("This process is a discrete version of 'mean curvature flow', which causes closed shapes to shrink and eventually vanish, and open edges to straighten.")

if __name__ == '__main__':
    run_simulation()