import numpy as np

def print_image(image):
    """Prints a binary image using '#' for 1 and '.' for 0."""
    for row in image:
        # The prompt asks to output each number, so we print the numeric representation
        # and a more visual character representation.
        numeric_row = "".join([str(p) for p in row])
        char_row = "".join(['#' if p == 1 else ' ' for p in row])
        print(f"{char_row}    ({numeric_row})")
    print("-" * image.shape[1])

def apply_local_median(image, radius):
    """
    Applies a local median filter to the image using a circular neighborhood.
    For each pixel, it calculates the median of the values in the neighborhood
    defined by the radius.
    """
    height, width = image.shape
    # Create a copy to store the new image state
    new_image = np.copy(image)
    
    # Create a grid of coordinates to define the circular neighborhood
    y, x = np.ogrid[-radius:radius+1, -radius:radius+1]
    mask = x**2 + y**2 <= radius**2

    for r in range(height):
        for c in range(width):
            # Define bounding box for the neighborhood to handle image borders
            r_min, r_max = max(0, r - radius), min(height, r + radius + 1)
            c_min, c_max = max(0, c - radius), min(width, c + radius + 1)

            # Extract the patch from the image and the corresponding mask part
            image_patch = image[r_min:r_max, c_min:c_max]
            
            mask_r_min = radius - (r - r_min)
            mask_r_max = radius + (r_max - r)
            mask_c_min = radius - (c - c_min)
            mask_c_max = radius + (c_max - c)
            
            final_mask = mask[mask_r_min:mask_r_max, mask_c_min:mask_c_max]
            
            # Get neighborhood values and calculate the median
            neighborhood_values = image_patch[final_mask]
            median_val = np.median(neighborhood_values)
            
            # For a binary image, if the median is 0.5 (equal 0s and 1s),
            # np.round will round to the nearest even integer, which is 0.
            # To be more deterministic, we can use int(median_val + 0.5) to always round up.
            # This makes the filter favor "filling" over "eroding" in ambiguous cases.
            new_image[r, c] = int(median_val + 0.5)
            
    return new_image

def run_simulation():
    """Sets up and runs the simulation."""
    # --- Simulation Parameters (Numbers from the 'equation') ---
    N = 30      # Image size: N x N
    delta = 3   # Radius of the median filter neighborhood
    timesteps_to_show = [0, 2, 5, 10, 15]

    print(f"Running simulation with parameters:")
    print(f"Image Size (N): {N}")
    print(f"Filter Radius (delta): {delta}")
    print("-" * N)

    # --- Initial Image ---
    # Create an initial image: a 'C' shape
    image = np.zeros((N, N), dtype=int)
    image[8:22, 8:22] = 1   # A solid square
    image[11:19, 11:19] = 0 # Hollow it out
    image[13:17, 15:22] = 0 # Add a gap to make it a 'C'

    current_image = np.copy(image)

    # --- Run Iterations ---
    max_iter = max(timesteps_to_show)
    for t in range(max_iter + 1):
        if t in timesteps_to_show:
            print(f"--- Image at t = {t} ---")
            print_image(current_image)
        
        # Stop if the image becomes completely uniform
        if np.all(current_image == 0) or np.all(current_image == 1):
            print(f"Image stabilized to a uniform state at t = {t}.")
            break
            
        current_image = apply_local_median(current_image, delta)

if __name__ == '__main__':
    run_simulation()