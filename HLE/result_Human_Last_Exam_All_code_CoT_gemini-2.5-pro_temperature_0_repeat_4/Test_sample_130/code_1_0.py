import numpy as np
from scipy.ndimage import median_filter
from skimage.morphology import disk

def print_image_ascii(img, downscale_factor=5):
    """
    Prints a downscaled ASCII art representation of the binary image.
    """
    h, w = img.shape
    # Ensure we don't have a zero step
    row_step = max(1, downscale_factor)
    col_step = max(1, downscale_factor * 2) # Adjust for character aspect ratio

    for r in range(0, h, row_step):
        row_str = ""
        for c in range(0, w, col_step):
            block = img[r:r+row_step, c:c+col_step]
            # Use mean to decide if the block is mostly 'on' or 'off'
            if np.mean(block) > 0.5:
                row_str += "##"
            else:
                row_str += "  "
        print(row_str)

def simulate_median_filter_evolution():
    """
    Simulates the evolution of a binary image under a repeated median filter
    and prints the state at different timesteps.
    """
    # --- Parameters ---
    # Image size (N x N)
    N = 100
    # Filter radius (delta)
    delta = 5
    # Timesteps to display
    timesteps_to_show = [0, 10, 30, 60]
    max_iterations = max(timesteps_to_show) + 1

    # --- Initialization ---
    # Create a black image
    image = np.zeros((N, N), dtype=np.uint8)
    # Add a white square in the middle
    start, end = N // 4, 3 * N // 4
    image[start:end, start:end] = 1

    # Create the circular neighborhood for the median filter
    # This corresponds to ||x-y|| <= delta
    footprint = disk(delta)

    print("--- Simulation of Median Filter Evolution ---")
    print(f"Image Size (N): {N}")
    print(f"Filter Radius (delta): {delta}")
    print("The operator is f(I) = Median(I) over a circular neighborhood.")
    print("We will observe the state of a white square at different timesteps (t).\n")

    # --- Simulation Loop ---
    for t in range(max_iterations):
        if t in timesteps_to_show:
            print("-" * 40)
            print(f"Timestep t = {t}")
            # The "equation" is the total number of white pixels
            total_white_pixels = np.sum(image)
            print(f"Total white pixels (area): {total_white_pixels}")
            print("Image shape:")
            print_image_ascii(image, downscale_factor=N//20)
            print("-" * 40 + "\n")

        # Apply the local median filter
        image = median_filter(image, footprint=footprint)

if __name__ == '__main__':
    simulate_median_filter_evolution()